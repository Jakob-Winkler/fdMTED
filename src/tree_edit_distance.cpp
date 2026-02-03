#include <cstdint>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <utility>

#include <Rcpp.h>

using namespace Rcpp;

// ChatGPT code
// ---- Small runtime-check helpers (fail fast instead of crashing) ----
static inline void check_index(std::size_t i, std::size_t n, const char* what) {
  if (i >= n) Rcpp::stop("%s: index %llu out of bounds (size %llu)",
                        what,
                        (unsigned long long)i,
                        (unsigned long long)n);
}
static inline void check_nonneg(int i, const char* what) {
  if (i < 0) Rcpp::stop("%s: negative index %d", what, i);
}
static inline void check_same_length(R_xlen_t a, R_xlen_t b, const char* what_a, const char* what_b) {
  if (a != b) Rcpp::stop("Length mismatch: %s has length %lld but %s has length %lld",
                        what_a, (long long)a, what_b, (long long)b);
}
static inline void validate_tree_vectors(const IntegerVector& parent,
                                        const IntegerVector& c1,
                                        const IntegerVector& c2,
                                        const NumericVector& w,
                                        const char* which) {
  check_same_length(parent.size(), c1.size(), "parent_vec", "child_1_vec");
  check_same_length(parent.size(), c2.size(), "parent_vec", "child_2_vec");
  check_same_length(parent.size(), w.size(),  "parent_vec", "weight_vec");

  const int n = parent.size();
  bool saw_root = false;

  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(parent[i])) saw_root = true;

    if (!IntegerVector::is_na(c1[i])) {
      if (c1[i] < 0 || c1[i] >= n) Rcpp::stop("%s: child_1_vec[%d]=%d out of range [0,%d)", which, i, c1[i], n);
    }
    if (!IntegerVector::is_na(c2[i])) {
      if (c2[i] < 0 || c2[i] >= n) Rcpp::stop("%s: child_2_vec[%d]=%d out of range [0,%d)", which, i, c2[i], n);
    }

    if (IntegerVector::is_na(c1[i]) != IntegerVector::is_na(c2[i])) {
      Rcpp::stop("%s: node %d has only one child set (child_1 NA? %d, child_2 NA? %d). Expected both NA (leaf) or both non-NA (internal).",
                 which, i,
                 (int)IntegerVector::is_na(c1[i]),
                 (int)IntegerVector::is_na(c2[i]));
    }
  }
  if (!saw_root) Rcpp::stop("%s: no root found (no NA in parent_vec)", which);
}

// Bit Vector OR (ChatGPT code)



namespace bitops {

// mask for bits [lo, hi) within a 64-bit word, where 0 <= lo < hi <= 64
static inline uint64_t mask_in_word(unsigned lo, unsigned hi) {
  uint64_t left  = (lo == 0)  ? ~0ULL : (~0ULL << lo);
  uint64_t right = (hi == 64) ? ~0ULL : ((1ULL << hi) - 1ULL);
  return left & right;
}

// set bits outside [0, len) to 0 in the last word of a packed array
static inline void mask_tail(uint64_t* a, std::size_t n_words, std::size_t len_bits) {
  if (n_words == 0) return;
  unsigned rem = static_cast<unsigned>(len_bits & 63ULL);
  if (rem == 0) return;
  uint64_t m = (1ULL << rem) - 1ULL;
  a[n_words - 1] &= m;
}

static inline void set_bit(uint64_t* words, std::size_t bit) {
  words[bit >> 6] |= (1ULL << (bit & 63));
}
static inline void clear_bit(uint64_t* words, std::size_t bit) {
  words[bit >> 6] &= ~(1ULL << (bit & 63));
}
static inline bool test_bit(const uint64_t* words, std::size_t bit) {
  return ((words[bit >> 6] >> (bit & 63)) & 1ULL) != 0ULL;
}

} // namespace bitops

class LeafConstraints {
public:
  LeafConstraints(std::size_t num_vecs, std::size_t num_bits)
    : V(num_vecs),
      L(num_bits),
      W((num_bits + 63) / 64),
      data(V * W, 0ULL),
      edit(W, 0ULL) {}

  std::size_t num_vecs() const { return V; }
  std::size_t num_bits() const { return L; }
  std::size_t words_per_vec() const { return W; }
  
  const std::vector<uint64_t>& edit_storage() const { return edit; }
  

  uint64_t* row(std::size_t vec_id) { return data.data() + vec_id * W; }
  const uint64_t* row(std::size_t vec_id) const { return data.data() + vec_id * W; }

  // --- single-bit ops on leaf_constraints ---
  void leaf_set(std::size_t vec_id, std::size_t bit) {
    check_index(vec_id, V, "LeafConstraints::leaf_set vec_id");
    check_index(bit,   L, "LeafConstraints::leaf_set bit");
    bitops::set_bit(row(vec_id), bit);
  }
  void leaf_clear(std::size_t vec_id, std::size_t bit) {
    check_index(vec_id, V, "LeafConstraints::leaf_clear vec_id");
    check_index(bit,   L, "LeafConstraints::leaf_clear bit");
    bitops::clear_bit(row(vec_id), bit);
  }
  bool leaf_test(std::size_t vec_id, std::size_t bit) const {
    check_index(vec_id, V, "LeafConstraints::leaf_test vec_id");
    check_index(bit,   L, "LeafConstraints::leaf_test bit");
    return bitops::test_bit(row(vec_id), bit);
  }

  // --- single-bit ops on edit_vec ---
  void edit_set(std::size_t bit) {
    check_index(bit, L, "LeafConstraints::edit_set bit");
    bitops::set_bit(edit.data(), bit);
  }
  void edit_clear(std::size_t bit) {
    check_index(bit, L, "LeafConstraints::edit_clear bit");
    bitops::clear_bit(edit.data(), bit);
  }
  bool edit_test(std::size_t bit) const {
    check_index(bit, L, "LeafConstraints::edit_test bit");
    return bitops::test_bit(edit.data(), bit);
  }

  // --- exact OR: leaf[vec][bit_lo,bit_hi) |= edit[bit_lo,bit_hi) ---
  void or_edit_into_leaf_range(std::size_t vec_id,
                              std::size_t bit_lo,
                              std::size_t bit_hi) {
    check_index(vec_id, V, "LeafConstraints::or_edit_into_leaf_range vec_id");
    if (bit_lo > L) bit_lo = L;
    if (bit_hi > L) bit_hi = L;
    if (bit_lo >= bit_hi) return;

    uint64_t* leaf = row(vec_id);
    const uint64_t* e = edit.data();

    const std::size_t w0 = bit_lo >> 6;
    const std::size_t w1 = (bit_hi - 1) >> 6;
    const unsigned b0 = static_cast<unsigned>(bit_lo & 63ULL);
    const unsigned b1 = static_cast<unsigned>(bit_hi & 63ULL);

    if (w0 == w1) {
      unsigned hi = (b1 == 0) ? 64U : b1;
      uint64_t m = bitops::mask_in_word(b0, hi);
      leaf[w0] |= (e[w0] & m);
      return;
    }

    // first partial word
    leaf[w0] |= (e[w0] & bitops::mask_in_word(b0, 64U));

    // full middle words
    for (std::size_t w = w0 + 1; w < w1; ++w) leaf[w] |= e[w];

    // last partial word
    unsigned hi = (b1 == 0) ? 64U : b1;
    leaf[w1] |= (e[w1] & bitops::mask_in_word(0U, hi));
  }

  // --- exact shift in edit_vec window ---
  void edit_shift_window(std::size_t bit_lo,
                         std::size_t bit_hi,
                         std::size_t shift_bits,
                         int dir) {
    if (bit_lo > L) bit_lo = L;
    if (bit_hi > L) bit_hi = L;
    if (bit_lo >= bit_hi) return;

    const std::size_t len = bit_hi - bit_lo;
    if (shift_bits >= len) {
      zero_window(bit_lo, bit_hi);
      return;
    }

    const std::size_t w_start = bit_lo >> 6;
    const std::size_t w_end   = (bit_hi - 1) >> 6;
    const std::size_t aw = w_end - w_start + 1;
    const unsigned offset = static_cast<unsigned>(bit_lo & 63ULL);

    std::vector<uint64_t> buf(aw, 0ULL);

    for (std::size_t i = 0; i < aw; ++i) {
      uint64_t a = edit[w_start + i];
      uint64_t b = (offset && w_start + i + 1 < W) ? edit[w_start + i + 1] : 0ULL;
      buf[i] = offset ? ((a >> offset) | (b << (64 - offset))) : a;
    }

    bitops::mask_tail(buf.data(), aw, len);

    if (dir > 0) shift_left_buf(buf.data(), aw, len, shift_bits);
    else         shift_right_buf(buf.data(), aw, len, shift_bits);

    write_back_window(buf.data(), aw, bit_lo, bit_hi);
  }

private:
  std::size_t V, L, W;
  std::vector<uint64_t> data;
  std::vector<uint64_t> edit;

  static void shift_left_buf(uint64_t* buf,
                             std::size_t n,
                             std::size_t len,
                             std::size_t s) {
    if (!s) return;
    std::size_t ws = s >> 6;
    unsigned bs = static_cast<unsigned>(s & 63ULL);

    if (ws) {
      for (std::size_t i = n; i-- > 0;) {
        buf[i] = (i >= ws) ? buf[i - ws] : 0ULL;
      }
    }
    if (bs) {
      for (std::size_t i = n; i-- > 0;) {
        uint64_t c = (i > 0) ? (buf[i - 1] >> (64 - bs)) : 0ULL;
        buf[i] = (buf[i] << bs) | c;
      }
    }
    bitops::mask_tail(buf, n, len);
  }

  static void shift_right_buf(uint64_t* buf,
                              std::size_t n,
                              std::size_t len,
                              std::size_t s) {
    if (!s) return;
    std::size_t ws = s >> 6;
    unsigned bs = static_cast<unsigned>(s & 63ULL);

    if (ws) {
      for (std::size_t i = 0; i < n; ++i) {
        buf[i] = (i + ws < n) ? buf[i + ws] : 0ULL;
      }
    }
    if (bs) {
      for (std::size_t i = 0; i < n; ++i) {
        uint64_t c = (i + 1 < n) ? (buf[i + 1] << (64 - bs)) : 0ULL;
        buf[i] = (buf[i] >> bs) | c;
      }
    }
    bitops::mask_tail(buf, n, len);
  }

  void zero_window(std::size_t bit_lo, std::size_t bit_hi) {
    uint64_t* e = edit.data();
    std::size_t w0 = bit_lo >> 6;
    std::size_t w1 = (bit_hi - 1) >> 6;
    unsigned b0 = static_cast<unsigned>(bit_lo & 63ULL);
    unsigned b1 = static_cast<unsigned>(bit_hi & 63ULL);

    if (w0 == w1) {
      unsigned hi = (b1 == 0) ? 64U : b1;
      e[w0] &= ~bitops::mask_in_word(b0, hi);
      return;
    }

    e[w0] &= ~bitops::mask_in_word(b0, 64U);
    for (std::size_t w = w0 + 1; w < w1; ++w) e[w] = 0ULL;
    unsigned hi = (b1 == 0) ? 64U : b1;
    e[w1] &= ~bitops::mask_in_word(0U, hi);
  }

  void write_back_window(const uint64_t* buf,
                         std::size_t aw,
                         std::size_t bit_lo,
                         std::size_t bit_hi) {
    uint64_t* e = edit.data();
    std::size_t w_start = bit_lo >> 6;
    std::size_t w_end   = (bit_hi - 1) >> 6;
    unsigned offset = static_cast<unsigned>(bit_lo & 63ULL);

    for (std::size_t w = w_start; w <= w_end; ++w) {
      std::size_t i = w - w_start;
      uint64_t nw;

      if (!offset) {
        nw = (i < aw) ? buf[i] : 0ULL;
      } else {
        uint64_t lo = (i < aw) ? (buf[i] << offset) : 0ULL;
        uint64_t ca = (i > 0) ? (buf[i - 1] >> (64 - offset)) : 0ULL;
        nw = lo | ca;
      }

      std::size_t wb0 = w * 64;
      unsigned lo = static_cast<unsigned>((bit_lo > wb0) ? (bit_lo - wb0) : 0);
      unsigned hi = static_cast<unsigned>((bit_hi < wb0 + 64) ? (bit_hi - wb0) : 64);

      uint64_t m = (lo == hi) ? 0ULL : bitops::mask_in_word(lo, hi);
      e[w] = (e[w] & ~m) | (nw & m);
    }
  }
};

// ChatGPT code
static inline void set_len_inplace(Rcpp::NumericVector& v, R_xlen_t n) {
  SEXP s = v;                 // same underlying object
  s = Rf_lengthgets(s, n);    // change *logical* length (no realloc if n <= TRUELENGTH)
  v = Rcpp::NumericVector(s); // rebind wrapper (safe)
}

// Transforms Tree into DFS representation 
static int initial_tree_parse(std::vector<std::pair<int, int>>& vistStack,
                              std::vector<int>& depthStack,
                              std::vector<int>& nodeStack,
                              std::vector<int>& tin,
                              std::vector<int>& tout,
                              IntegerVector& parent_vec,
                              IntegerVector& child_1_vec,
                              IntegerVector& child_2_vec) {
  int timer = 0;
  int hight = 0;
  int counter = 0;
  int n_leafs = 0;

  std::vector<int> counterStack;
  counterStack.reserve((std::size_t)child_1_vec.size());

  // find root
  for (int i = 0; i < parent_vec.size(); ++i) {
    if (IntegerVector::is_na(parent_vec[i])) {
      vistStack.push_back({i, 0});
      break;
    }
  }

  while (!vistStack.empty()) {
    int v = vistStack.back().first;
    if (v < 0 || v >= parent_vec.size()) Rcpp::stop("initial_tree_parse: v=%d out of range", v);
    int& state = vistStack.back().second;

    if (state == 0) {
      // enter
      depthStack.push_back(hight);
      nodeStack.push_back(v); // Thats bad
      tin.push_back(timer++);

      // leaf
      if (IntegerVector::is_na(child_1_vec[v])) {
        tout[counter++] = timer - 1;
        vistStack.pop_back();
        hight--;
        n_leafs++;
        continue;
      }

      counterStack.push_back(counter++);
      state = 1;
      hight++;
      vistStack.push_back({child_1_vec[v], 0});
      continue;
    }

    if (state == 1) {
      state = 2;
      if (IntegerVector::is_na(child_2_vec[v])) Rcpp::stop("initial_tree_parse: internal node %d has NA child_2", v);
      int c2 = child_2_vec[v];
      if (c2 < 0 || c2 >= parent_vec.size()) Rcpp::stop("initial_tree_parse: child_2_vec[%d]=%d out of range", v, c2);
      hight++;
      vistStack.push_back({c2, 0});
      continue;
    }

    // exit
    hight--;
    tout[counterStack.back()] = timer - 1;
    counterStack.pop_back();
    vistStack.pop_back();
  }

  return n_leafs;
}

// Tree data structure (is all the information nessary to run solver)
struct tree_BLP_optim {
  LeafConstraints lc{0, 0};
  std::vector<double> subtrees;
  std::vector<double> path_result;
  std::vector<int> depth;
  std::vector<int> tin;
  std::vector<int> tout;
  int n_nodes = 0;
  int n_paths = 0;
  int n_leafs = 0;
};

// ChatGPT code
// ---- Expand one packed bitvector (W words) to 0/1 ints (length L) ----
static inline IntegerVector unpack_bits_to_int01(const uint64_t* words, std::size_t L) {
  IntegerVector out(L);
  std::size_t full_words = L >> 6;         // L / 64
  unsigned rem_bits = (unsigned)(L & 63);  // L % 64
  
  std::size_t idx = 0;
  for (std::size_t w = 0; w < full_words; ++w) {
    uint64_t x = words[w];
    for (unsigned b = 0; b < 64; ++b, ++idx) {
      out[idx] = (int)((x >> b) & 1ULL);
    }
  }
  if (rem_bits) {
    uint64_t x = words[full_words];
    for (unsigned b = 0; b < rem_bits; ++b, ++idx) {
      out[idx] = (int)((x >> b) & 1ULL);
    }
  }
  return out;
}


// ---- Expand leaf_constraints (V x L) to IntegerMatrix of 0/1 ----
static inline IntegerMatrix unpack_leaf_to_matrix01(const LeafConstraints& lc) {
  const std::size_t V = lc.num_vecs();
  const std::size_t L = lc.num_bits();
  
  IntegerMatrix M((int)V, (int)L);
  
  for (std::size_t i = 0; i < V; ++i) {
    const uint64_t* row_words = lc.row(i);
    
    std::size_t full_words = L >> 6;
    unsigned rem_bits = (unsigned)(L & 63);
    
    std::size_t col = 0;
    
    for (std::size_t w = 0; w < full_words; ++w) {
      uint64_t x = row_words[w];
      for (unsigned b = 0; b < 64; ++b, ++col) {
        M((int)i, (int)col) = (int)((x >> b) & 1ULL);
      }
    }
    if (rem_bits) {
      uint64_t x = row_words[full_words];
      for (unsigned b = 0; b < rem_bits; ++b, ++col) {
        M((int)i, (int)col) = (int)((x >> b) & 1ULL);
      }
    }
  }
  
  return M;
}

// ---- Convert tree_BLP_optim -> R list with 0/1 bit vectors ----
static inline List tree_BLP_optim_to_list_bits01(const tree_BLP_optim& out) {
  const std::size_t L = out.lc.num_bits();
  
  IntegerMatrix leaf01 = unpack_leaf_to_matrix01(out.lc);
  
  // edit_vec: expand from its packed words
  const std::vector<uint64_t>& edit_words = out.lc.edit_storage();
  IntegerVector edit01 = unpack_bits_to_int01(edit_words.data(), L);
  
  return List::create(
    _["leaf_constraints"] = leaf01,     // V x L matrix of 0/1
    _["edit_vec"] = edit01,             // length L vector of 0/1
    _["subtrees"] = wrap(out.subtrees),
    _["path_result"] = wrap(out.path_result),
    _["depth"] = wrap(out.depth),
    _["tin"] = wrap(out.tin),
    _["tout"] = wrap(out.tout)
  );
}
// Transforms input tree to tree_BLP_optim representation
static tree_BLP_optim tree_preprosses(IntegerVector parent_vec,
                                     IntegerVector child_1_vec,
                                     IntegerVector child_2_vec,
                                     NumericVector weight_vec) {
  validate_tree_vectors(parent_vec, child_1_vec, child_2_vec, weight_vec, "tree_preprosses");
  
  // Define and reserve space for the output
  const R_xlen_t n_nodes = child_1_vec.size();

  tree_BLP_optim tree;
  tree.n_nodes = (int)n_nodes;

  std::vector<std::pair<int, int>> vistStack;
  vistStack.reserve((std::size_t)n_nodes);

  tree.depth.reserve((std::size_t)n_nodes);

  std::vector<int> nodeStack;
  nodeStack.reserve((std::size_t)n_nodes);

  tree.tin.reserve((std::size_t)n_nodes);
  tree.tout.resize((std::size_t)n_nodes);

  // DFS 
  int n_leafs = initial_tree_parse(vistStack, tree.depth, nodeStack, tree.tin, tree.tout,
                                  parent_vec, child_1_vec, child_2_vec);
  tree.n_leafs = n_leafs;
  const int n = (int)nodeStack.size();
  NumericVector weight_vec_reorder(n);
  for (int i = 0; i < n; ++i) {
    weight_vec_reorder[i] = weight_vec[nodeStack[i]];
  }

  int n_paths = 0;
  for (int i = 1; i < n_nodes; ++i) {
    n_paths += tree.depth[i];
  }

  tree.lc = LeafConstraints((std::size_t)n_leafs, (std::size_t)n_paths);
  tree.n_paths = n_paths;

  tree.path_result.reserve((std::size_t)n_paths);
  tree.subtrees.resize((std::size_t)n_nodes);

  int max_depth = *std::max_element(tree.depth.begin(), tree.depth.end());
  std::vector<double> path_stack;
  path_stack.reserve((std::size_t)max_depth);

  std::vector<int> subtree_start_stack;
  subtree_start_stack.reserve((std::size_t)max_depth);

  std::vector<int> subtree_bit_start((std::size_t)n_nodes);

  int leaf_start_bit;
  int leaf_end_bit;
  int cur_leaf = 0;
  int subtree_leaf_start;
  double subtree_weight;
  double tmp_stack_weight;

  // Here all the cost and constrains are calculated. Later they are used to with the result of the other tree to get the final cost_vec and constrains matrix
  // This is in DFS order
  for (int i = 1; i < n_nodes; ++i) {
    path_stack.push_back(weight_vec_reorder[i]);

    subtree_bit_start[i] = subtree_bit_start[i - 1] + tree.depth[i - 1];

    if (subtree_bit_start[i] < 0 || subtree_bit_start[i] >= tree.n_paths)
      Rcpp::stop("tree_preprosses: subtree_bit_start[%d]=%d out of range [0,%d)", i, subtree_bit_start[i], tree.n_paths);
    tree.lc.edit_set((std::size_t)subtree_bit_start[i]);

    tmp_stack_weight = 0;
    
    for (int path_stack_idx = path_stack.size()-1; path_stack_idx>=0;path_stack_idx--) {
      tmp_stack_weight += path_stack[path_stack_idx];
      tree.path_result.push_back(tmp_stack_weight);
      
    }
    
    // We go until there is a leaf
    if (tree.tin[i] == tree.tout[i]) {
      subtree_weight = weight_vec_reorder[i];

      leaf_start_bit = subtree_bit_start[i];
      leaf_end_bit = (int)tree.path_result.size();

      
      tree.lc.leaf_set((std::size_t)cur_leaf, (std::size_t)leaf_start_bit);
      tree.lc.edit_set((std::size_t)leaf_start_bit);
      
      tree.lc.edit_shift_window((std::size_t)leaf_start_bit,
                                (std::size_t)leaf_end_bit,
                                1,
                                1);

      subtree_leaf_start = cur_leaf;

      // From the leaf we fall back up the subtree
      for (int j = i - 1; j >= 0; --j) { 
        
        // We are on the right side of the subtree, so we "merge"
        if (tree.tout[i] == tree.tout[j]) {
          
          tree.subtrees[j] += subtree_weight;
          subtree_weight = tree.subtrees[j] + weight_vec_reorder[j];
          path_stack.pop_back();

          if(j == 0) break; // on fall back, I the root is reached stop 
          leaf_start_bit = subtree_bit_start[j];
          

          subtree_leaf_start = subtree_start_stack.back();
          subtree_start_stack.pop_back();
          
          for (int leaf = subtree_leaf_start; leaf <= cur_leaf; ++leaf) {
            tree.lc.or_edit_into_leaf_range((std::size_t)leaf,
                                           (std::size_t)leaf_start_bit,
                                           (std::size_t)leaf_end_bit);
          }

          tree.lc.edit_shift_window((std::size_t)leaf_start_bit,
                                   (std::size_t)leaf_end_bit,
                                   1,
                                   1);
        }

        // We are on the left side of a subtree, so we safe the possiton for the right pass
        if (tree.tout[i] < tree.tout[j]) {
          
          tree.subtrees[j] = subtree_weight;
          path_stack.pop_back();

          subtree_start_stack.push_back(subtree_leaf_start);
          
          break;
        }
      }

      cur_leaf++;
    }
  }

  
  return tree;
}

// For each subtree gets the relevant nodes/leafs/paths
static void get_sel_vec(const tree_BLP_optim& tree,
                        int idx_tree,
                        std::vector<int>& tree_path_sel,
                        std::vector<int>& tree_path_node_sel,
                        std::pair<int, int>& tree_leaf_sel) {
  check_nonneg(idx_tree, "get_sel_vec idx_tree");
  check_index((std::size_t)idx_tree, tree.depth.size(), "get_sel_vec idx_tree");
  int offset = 0;
  int leafs_seen = 0;
  std::size_t i = 0;

  for (; i < (std::size_t)idx_tree; ++i) {
    offset += tree.depth[i];
    if (tree.tin[i] == tree.tout[i]) leafs_seen++;
  }
  tree_leaf_sel.first = leafs_seen;

  for (; i < tree.depth.size(); ++i) {
    if (tree.tout[i] > tree.tout[(std::size_t)idx_tree]) break;

    if (tree.tin[i] == tree.tout[i]) leafs_seen++;

    for (int j = 0; j < (tree.depth[i] - tree.depth[(std::size_t)idx_tree]); ++j) { // For subTree of node A, node A is not part of edit set
      tree_path_sel.push_back(offset + j);
      tree_path_node_sel.push_back((int)i);
    }
    offset += tree.depth[i];
  }

  tree_leaf_sel.second = leafs_seen;
}

// Using the selected var from get_sel_vec construct cost and constrains to be solved
static void gen_BLP(const tree_BLP_optim& tree1,
                    const tree_BLP_optim& tree2,
                    const std::vector<int>& tree1_path_sel,
                    const std::vector<int>& tree1_path_node_sel,
                    const std::pair<int, int>& tree1_leaf_range,
                    const std::vector<int>& tree2_path_sel,
                    const std::vector<int>& tree2_path_node_sel,
                    const std::pair<int, int>& tree2_leaf_range,
                    NumericVector obj_vec,
                    NumericVector const_Matrix,
                    const int& const_matrix_ncol,
                    const NumericMatrix& W) {
  
  const std::size_t n1 = tree1_path_sel.size();
  const std::size_t n2 = tree2_path_sel.size();
  
  const int n_leaf1 = tree1_leaf_range.second - tree1_leaf_range.first;
  const int n_leaf2 = tree2_leaf_range.second - tree2_leaf_range.first;
  const int n_rows  = n_leaf1 + n_leaf2;

  
  // Keep ORIGINAL indexing convention:
  // const_Matrix[row * const_matrix_ncol + col]
  auto cm_idx = [n_rows](int row, int col) -> R_xlen_t {
    return static_cast<R_xlen_t>(row) + static_cast<R_xlen_t>(n_rows) * static_cast<R_xlen_t>(col);
  };
  
  for (std::size_t i = 0; i < n1; ++i) {
    const std::size_t t1_path_idx = static_cast<std::size_t>(tree1_path_sel[i]);
    const std::size_t t1_node_idx = static_cast<std::size_t>(tree1_path_node_sel[i]);
    
    check_index(t1_path_idx, tree1.path_result.size(), "gen_BLP tree1_path_sel");
    check_index(t1_node_idx, tree1.subtrees.size(),    "gen_BLP tree1_path_node_sel");
    
    const double t1_path_val = tree1.path_result[t1_path_idx];
    const double t1_subtree  = tree1.subtrees[t1_node_idx];
    
    for (std::size_t j = 0; j < n2; ++j) {
      const std::size_t t2_path_idx = static_cast<std::size_t>(tree2_path_sel[j]);
      const std::size_t t2_node_idx = static_cast<std::size_t>(tree2_path_node_sel[j]);
      
      check_index(t2_path_idx, tree2.path_result.size(), "gen_BLP tree2_path_sel");
      check_index(t2_node_idx, tree2.subtrees.size(),    "gen_BLP tree2_path_node_sel");
      
      const double t2_path_val = tree2.path_result[t2_path_idx];
      const double t2_subtree  = tree2.subtrees[t2_node_idx];
      
      const int col = static_cast<int>(i * n2 + j);
      
      // Objective (unchanged)
      obj_vec[col] =
        std::abs(t1_path_val - t2_path_val)
        - t1_path_val - t2_path_val
      - t1_subtree - t2_subtree
      + W(tree1_path_node_sel[i], tree2_path_node_sel[j]);
      
      // Constraints: first tree1 leaves, then tree2 leaves (same as original)
      for (int r = 0; r < n_leaf1; ++r) {
        const std::size_t leaf = static_cast<std::size_t>(tree1_leaf_range.first + r);
        const_Matrix[cm_idx(r, col)] =
          static_cast<int>(tree1.lc.leaf_test(leaf, t1_path_idx));
      }
      
      for (int r = 0; r < n_leaf2; ++r) {
        const std::size_t leaf = static_cast<std::size_t>(tree2_leaf_range.first + r);
        const_Matrix[cm_idx(n_leaf1 + r, col)] =
          static_cast<int>(tree2.lc.leaf_test(leaf, t2_path_idx));
      }
    }
  }
  
}


// Takes subtrees and returns the distance (it calles the solver function)
static double calc_BLP(const tree_BLP_optim& tree1,
                       const tree_BLP_optim& tree2,
                       int idx_tree1,
                       int idx_tree2,
                       const NumericMatrix& W,
                       Function BLP_solver_call,
                       NumericVector& obj_buf,
                       NumericVector& cm_buf) {

  std::vector<int> tree1_path_sel;
  tree1_path_sel.reserve(tree1.path_result.size());
  std::vector<int> tree1_path_node_sel;
  tree1_path_node_sel.reserve(tree1.path_result.size());
  std::pair<int, int> tree1_leaf_range;

  get_sel_vec(tree1, idx_tree1, tree1_path_sel, tree1_path_node_sel, tree1_leaf_range);

  std::vector<int> tree2_path_sel;
  tree2_path_sel.reserve(tree2.path_result.size());
  std::vector<int> tree2_path_node_sel;
  tree2_path_node_sel.reserve(tree2.path_result.size());
  std::pair<int, int> tree2_leaf_range;

  get_sel_vec(tree2, idx_tree2, tree2_path_sel, tree2_path_node_sel, tree2_leaf_range);

  double BLP_offset = tree1.subtrees[(std::size_t)idx_tree1] + tree2.subtrees[(std::size_t)idx_tree2];

  int n_cols = (int)(tree1_path_sel.size() * tree2_path_sel.size());
  if (n_cols <= 0) Rcpp::stop("calc_BLP: empty selection");
  
  const int n_leaf1 = tree1_leaf_range.second - tree1_leaf_range.first;
  const int n_leaf2 = tree2_leaf_range.second - tree2_leaf_range.first;
  const int n_rows  = n_leaf1 + n_leaf2;
  
  const R_xlen_t need_obj = (R_xlen_t)n_cols;
  const R_xlen_t need_cm  = (R_xlen_t)n_rows * (R_xlen_t)n_cols;

  
  set_len_inplace(obj_buf, need_obj);
  set_len_inplace(cm_buf,  need_cm);
  cm_buf.attr("dim") = Rcpp::Dimension(n_rows, n_cols);
 

  gen_BLP(tree1, tree2,
          tree1_path_sel, tree1_path_node_sel, tree1_leaf_range,
          tree2_path_sel, tree2_path_node_sel, tree2_leaf_range,
          obj_buf, cm_buf, n_cols, W); // ob_vec size is never set

  SEXP res = BLP_solver_call(obj_buf, cm_buf, BLP_offset, idx_tree1, idx_tree2, W);
  if (!Rf_isReal(res)) Rcpp::stop("BLP_solver must return a numeric scalar");
  NumericVector out(res);
  if (out.size() < 1) Rcpp::stop("BLP_solver returned empty result");
  return out[0];
}

// Here the Bottom up algorithm is done
// The order is by subtree
// The root distance is returned
static NumericVector edit_dist_subtree_loop(const tree_BLP_optim& tree_1,
                                           const tree_BLP_optim& tree_2) {
  Function BLP_solver_call("BLP_solver_Cpp");

  int nrow = (int)tree_1.depth.size();
  int ncol = (int)tree_2.depth.size();
  if (nrow <= 0 || ncol <= 0) Rcpp::stop("edit_dist_subtree_loop: empty tree(s)");
  NumericMatrix W(nrow, ncol);

  // ---- allocate max buffers ONCE ----
  const R_xlen_t max_cols =
    (R_xlen_t)tree_1.path_result.size() * (R_xlen_t)tree_2.path_result.size();
  
  const R_xlen_t max_rows =
    (R_xlen_t)tree_1.lc.num_vecs() + (R_xlen_t)tree_2.lc.num_vecs();
  
  NumericVector obj_buf(Rcpp::no_init(max_cols));
  NumericVector cm_buf(Rcpp::no_init(max_cols * max_rows));
  
  // Tree 1 loop
  for (std::size_t idx_tree1 = 0; idx_tree1 < tree_1.depth.size(); ++idx_tree1) {
    if (tree_1.tin[idx_tree1] == tree_1.tout[idx_tree1]) {
      for (int sub_idx_tree1 = (int)idx_tree1; 0 <= sub_idx_tree1; --sub_idx_tree1) {
        if (tree_1.tout[idx_tree1] == tree_1.tout[(std::size_t)sub_idx_tree1]) {
          // Tree 2 loop
          for (std::size_t idx_tree2 = 0; idx_tree2 < tree_2.depth.size(); ++idx_tree2) {
            if (tree_2.tin[idx_tree2] == tree_2.tout[idx_tree2]) {
              for (int sub_idx_tree2 = (int)idx_tree2; 0 <= sub_idx_tree2; --sub_idx_tree2) {
                if (tree_2.tout[idx_tree2] == tree_2.tout[(std::size_t)sub_idx_tree2]) {

                  // base cases
                  if (tree_1.tin[(std::size_t)sub_idx_tree1] == tree_1.tout[(std::size_t)sub_idx_tree1]) {
                    if (tree_2.tin[(std::size_t)sub_idx_tree2] == tree_2.tout[(std::size_t)sub_idx_tree2]) {
                      
                      continue;
                    } else {
                      W(sub_idx_tree1, sub_idx_tree2) = tree_2.subtrees[(std::size_t)sub_idx_tree2];
                      continue;
                    }
                  } else if (tree_2.tin[(std::size_t)sub_idx_tree2] == tree_2.tout[(std::size_t)sub_idx_tree2]) {
                    W(sub_idx_tree1, sub_idx_tree2) = tree_1.subtrees[(std::size_t)sub_idx_tree1];
                    continue;
                  } 

                  W(sub_idx_tree1, sub_idx_tree2) =
                    calc_BLP(tree_1, tree_2, sub_idx_tree1, sub_idx_tree2, W, BLP_solver_call, obj_buf, cm_buf);
                }
              }
            }
          }
        }
      }
    }
  }

  return NumericVector::create(W(0, 0));
}

// Function to be called in R
// [[Rcpp::export]]
NumericVector tree_edit_distance(IntegerVector parent_vec_1,
                                 IntegerVector child_1_vec_1,
                                 IntegerVector child_2_vec_1,
                                 NumericVector weight_vec_1,
                                 IntegerVector parent_vec_2,
                                 IntegerVector child_1_vec_2,
                                 IntegerVector child_2_vec_2,
                                 NumericVector weight_vec_2) {
  validate_tree_vectors(parent_vec_1, child_1_vec_1, child_2_vec_1, weight_vec_1, "tree_1");
  validate_tree_vectors(parent_vec_2, child_1_vec_2, child_2_vec_2, weight_vec_2, "tree_2");
  tree_BLP_optim tree_1 = tree_preprosses(parent_vec_1, child_1_vec_1, child_2_vec_1, weight_vec_1);
  tree_BLP_optim tree_2 = tree_preprosses(parent_vec_2, child_1_vec_2, child_2_vec_2, weight_vec_2);
  return edit_dist_subtree_loop(tree_1, tree_2);
}

// Helper funtion for debugging
// [[Rcpp::export]]
List get_prepross_tree(IntegerVector parent_vec_1,
                       IntegerVector child_1_vec_1,
                       IntegerVector child_2_vec_1,
                       NumericVector weight_vec_1) {
  tree_BLP_optim tree_1 = tree_preprosses(parent_vec_1,
                                          child_1_vec_1,
                                          child_2_vec_1,
                                          weight_vec_1);
  List out = tree_BLP_optim_to_list_bits01(tree_1);
  return out;
}
