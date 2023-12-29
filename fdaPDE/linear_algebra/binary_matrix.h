// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __BINARY_MATRIX_H__
#define __BINARY_MATRIX_H__

#include "../utils/symbols.h"
#include "../utils/assert.h"

#include <bitset>

namespace fdapde {
namespace core {

// forward declarations
template <int Rows, int Cols, typename E> class BinMtxBase;

// a space-time efficient implementation for a dense-storage matrix type holding binary values
template <int Rows, int Cols = Rows> class BinaryMatrix : public BinMtxBase<Rows, Cols, BinaryMatrix<Rows, Cols>> {
   public:
    using This = BinaryMatrix<Rows, Cols>;
    using Base = BinMtxBase<Rows, Cols, This>;
    using BitPackType = std::uintmax_t;
    static constexpr std::size_t PackSize = sizeof(std::uintmax_t) * 8; // number of bits in a packet
    static constexpr int ct_rows = Rows;
    static constexpr int ct_cols = Cols;
  
    template <typename T> struct is_dynamic_sized {
        static constexpr bool value = (T::ct_rows == Dynamic || T::ct_cols == Dynamic);
    };
    using StorageType = typename std::conditional<
      is_dynamic_sized<This>::value, std::vector<BitPackType>,
      std::array<BitPackType, 1 + (Rows * Cols) / PackSize>>::type;
    using Base::n_cols_;
    using Base::n_rows_;  
    // static default constructor
    fdapde_enable_constructor_if(!is_dynamic_sized, This) BinaryMatrix() : Base(Rows, Cols) {
        std::fill_n(data_.begin(), data_.size(), 0);
    };
    // dynamic sized constructors
    fdapde_enable_constructor_if(is_dynamic_sized, This) BinaryMatrix() : Base() {};
    fdapde_enable_constructor_if(is_dynamic_sized, This) BinaryMatrix(std::size_t n_rows, std::size_t n_cols) :
        Base(n_rows, n_cols) {
        data_.resize(1 + std::ceil((n_rows_ * n_cols_) / PackSize), 0);
    }
    // vector constructor
    fdapde_enable_constructor_if(is_dynamic_sized, This) BinaryMatrix(std::size_t n_rows) : BinaryMatrix(n_rows, 1) {
        fdapde_static_assert(Rows == Dynamic && Cols == 1, THIS_METHOD_IS_ONLY_FOR_VECTORS);
    }

    template <int Rows_ = Rows, int Cols_ = Cols>
    typename std::enable_if<Rows_ == Dynamic || Cols_ == Dynamic, void>::type
    resize(std::size_t rows, std::size_t cols) {
        n_rows_ = rows;
        n_cols_ = cols;
        data_ = std::vector<BitPackType>(1 + std::ceil((n_rows_ * n_cols_) / PackSize), 0);
    }

    // getters
    bool operator()(std::size_t i, std::size_t j) const {   // access to (i,j)-th element
        return (data_[pack_of(i, j)] & BitPackType(1) << ((i * Base::n_cols_ + j) % PackSize)) != 0;
    }
    template <int Rows_ = Rows, int Cols_ = Cols>   // vector-like (subscript) access
    typename std::enable_if<Rows_ == Dynamic && Cols_ == 1, bool>::type operator[](std::size_t i) const {
        return operator()(i, 0);
    }
    BitPackType bitpack(std::size_t i) const { return data_[i]; }
    // setters
    void set(std::size_t i, std::size_t j) {   // set (i,j)-th bit
        fdapde_assert(i < n_rows_ && j < n_cols_);
        data_[pack_of(i, j)] |= (BitPackType(1) << ((i * Base::n_cols_ + j) % PackSize));
    }
    void set(std::size_t i) {
        fdapde_static_assert(Rows == Dynamic && Cols == 1, THIS_METHOD_IS_ONLY_FOR_VECTORS);
        fdapde_assert(i < n_rows_);
	set(i, 0);
    }  
    void clear(std::size_t i, std::size_t j) {   // clear (i,j)-th bit (sets to 0)
        fdapde_assert(i < n_rows_ && j < n_cols_);
        data_[pack_of(i, j)] &= ~(BitPackType(1) << ((i * Base::n_cols_ + j) % PackSize));
    }
    void clear(std::size_t i) {
        fdapde_static_assert(Rows == Dynamic && Cols == 1, THIS_METHOD_IS_ONLY_FOR_VECTORS);
        fdapde_assert(i < n_rows_);
	clear(i, 0);
    }  
    // fast bitwise assignment from binary expression
    template <int Rows_, int Cols_, typename Rhs_> BinaryMatrix& operator=(const BinMtxBase<Rows_, Cols_, Rhs_>& rhs) {
        // !is_dynamic_sized \implies (Rows == Rows_ && Cols == Cols_)
        static_assert(is_dynamic_sized<This>::value || (Rows == Rows_ && Cols == Cols_));
	// perform assignment
        n_rows_ = rhs.rows();
        n_cols_ = rhs.cols();
        Base::n_bitpacks_ = rhs.bitpacks();
        if constexpr (is_dynamic_sized<This>::value) data_.resize(Base::n_bitpacks_, 0);
        for (std::size_t i = 0; i < rhs.bitpacks(); ++i) { data_[i] = rhs.bitpack(i); }
        return *this;
    }
   private:
    StorageType data_;
    // recover the byte-pack for the (i,j)-th element
    inline std::size_t pack_of(std::size_t i, std::size_t j) const { return (i * Base::n_cols_ + j) / PackSize; }
};

// unary biwise NOT of binary expression
template <int Rows, int Cols, typename XprTypeNested>
struct BinMtxNegationOp : public BinMtxBase<Rows, Cols, BinMtxNegationOp<Rows, Cols, XprTypeNested>> {
    using XprType = BinMtxNegationOp<Rows, Cols, XprTypeNested>;
    using Base = BinMtxBase<Rows, Cols, XprType>;
    using BitPackType = typename Base::BitPackType;
    const XprTypeNested& op_;

    BinMtxNegationOp(const XprTypeNested& op, int rows, int cols) : Base(rows, cols), op_(op) {};
    bool operator()(std::size_t i, std::size_t j) const { return !op_(i, j); }
    BitPackType bitpack(std::size_t i) const { return ~op_.bitpack(i); }   // bitwise not
};

// binary bitwise operations between binary expression
template <int Rows, int Cols, typename Lhs, typename Rhs, typename BinaryOperation>
struct BinMtxBinOp : public BinMtxBase<Rows, Cols, BinMtxBinOp<Rows, Cols, Lhs, Rhs, BinaryOperation>> {
    using XprType = BinMtxBinOp<Rows, Cols, Lhs, Rhs, BinaryOperation>;
    using Base = BinMtxBase<Rows, Cols, XprType>;
    using BitPackType = typename Base::BitPackType;
    const Lhs& op1_;
    const Rhs& op2_;
    BinaryOperation f_;

    BinMtxBinOp(const Lhs& op1, const Rhs& op2, BinaryOperation f) :
        Base(op1.rows(), op1.cols()), op1_(op1), op2_(op2), f_(f) {
        fdapde_assert(op1_.rows() == op2_.rows() && op1_.cols() == op2_.cols());
    }
    bool operator()(std::size_t i, std::size_t j) const { return f_(op1_(i, j), op2_(i, j)); }
    BitPackType bitpack(std::size_t i) const { return f_(op1_.bitpack(i), op2_.bitpack(i)); }
};

#define DEFINE_BITWISE_BINMTX_OPERATOR(OPERATOR, FUNCTOR)                                                              \
    template <int Rows, int Cols, typename OP1, typename OP2>                                                          \
    BinMtxBinOp<Rows, Cols, OP1, OP2, FUNCTOR> OPERATOR(                                                               \
      const BinMtxBase<Rows, Cols, OP1>& op1, const BinMtxBase<Rows, Cols, OP2>& op2) {                                \
        return BinMtxBinOp<Rows, Cols, OP1, OP2, FUNCTOR> {op1.get(), op2.get(), FUNCTOR()};                           \
    }

DEFINE_BITWISE_BINMTX_OPERATOR(operator&, std::bit_and<>)   // m1 & m2
DEFINE_BITWISE_BINMTX_OPERATOR(operator|, std::bit_or<> )   // m1 | m2
DEFINE_BITWISE_BINMTX_OPERATOR(operator^, std::bit_xor<>)   // m1 ^ m2

// dense-block of binary matrix
template <int BlockRows, int BlockCols, typename XprTypeNested>
class BinMtxBlock : public BinMtxBase<BlockRows, BlockCols, BinMtxBlock<BlockRows, BlockCols, XprTypeNested>> {
   private:
    using XprType = BinMtxBlock<BlockRows, BlockCols, XprTypeNested>;
    using Base = BinMtxBase<BlockRows, BlockCols, XprType>;
    using BitPackType = typename Base::BitPackType;
    static constexpr std::size_t PackSize = Base::PackSize;   // number of bits in a packet
    using Base::n_cols_;
    using Base::n_rows_;
    // internal data
    std::size_t start_row_, start_col_;
    XprTypeNested& xpr_;
   public:
    // row/column constructor
    BinMtxBlock(XprTypeNested& xpr, std::size_t i) :
        Base(BlockRows == 1 ? 1 : xpr.rows(), BlockCols == 1 ? 1 : xpr.cols()),
        xpr_(xpr),
        start_row_(BlockRows == 1 ? i : 0),
        start_col_(BlockCols == 1 ? i : 0) {
        fdapde_static_assert(BlockRows == 1 || BlockCols == 1, THIS_METHOD_IS_ONLY_FOR_ROW_AND_COLUMN_BLOCKS);
        fdapde_assert(i >= 0 && ((BlockRows == 1 && i < xpr_.rows()) || (BlockCols == 1 && i < xpr_.cols())));
    }
    // fixed-sized constructor
    BinMtxBlock(XprTypeNested& xpr, std::size_t start_row, std::size_t start_col) :
        Base(BlockRows, BlockCols), xpr_(xpr), start_row_(start_row), start_col_(start_col) {
        fdapde_static_assert(BlockRows != Dynamic && BlockCols != Dynamic, THIS_METHOD_IS_ONLY_FOR_FIXED_SIZE);
        fdapde_assert(
          start_row_ >= 0 && BlockRows >= 0 && start_row_ + BlockRows <= xpr_.rows() && start_col_ >= 0 &&
          BlockCols >= 0 && start_col_ + BlockCols <= xpr_.cols());
    }
    // dynamic-sized constructor
    BinMtxBlock(
      XprTypeNested& xpr, std::size_t start_row, std::size_t start_col, std::size_t block_rows,
      std::size_t block_cols) :
        Base(block_rows, block_cols), xpr_(xpr), start_row_(start_row), start_col_(start_col) {
        fdapde_assert(BlockRows == Dynamic || BlockCols == Dynamic);
        fdapde_assert(
          start_row_ >= 0 && start_row_ + block_rows <= xpr_.rows() && start_col_ >= 0 &&
          start_col_ + block_cols <= xpr_.cols());
    }

    bool operator()(std::size_t i, std::size_t j) const { return xpr_(i + start_row_, j + start_col_); }
    void set(std::size_t i, std::size_t j) { xpr_.set(i + start_row_, j + start_col_); }
    void clear(std::size_t i, std::size_t j) { xpr_.clear(i + start_row_, j + start_col_); }
    BitPackType bitpack(std::size_t i) const {
        BitPackType in = xpr_.bitpack(i);   // evaluate bitpack in a temporary
        BitPackType out = 0x0;
        // bitpack column offset
        std::size_t offset =
          start_row_ * xpr_.cols() + start_col_ * xpr_.rows() +        // bitpack offset
          (i == 0 ? 0 : xpr_.cols() - (i * PackSize % xpr_.cols()));   // possible bitpack misalignment offset
        std::size_t block_stride = xpr_.cols() - n_cols_;
	std::size_t base_offset = start_row_ * xpr_.cols();
        // extract bits from in packet according to block bit pattern
        BitPackType mask = 0x1;
        for (std::size_t i = 0; i < n_rows_; i++) {
            for (std::size_t j = 0; j < n_cols_; j++) {
                out |= (in & (mask << (offset + i * xpr_.cols() + j))) >> (base_offset + i * block_stride);
            }
        }
        return out;
    }
};

// visitors
template <typename XprType, typename Visitor> struct bitpack_visit {
    static constexpr std::size_t PackSize = XprType::PackSize;
    // apply visitor by cycling over bitpacks
    static inline void run(const XprType& xpr, Visitor& visitor) {
        std::size_t size = xpr.size();   // overall number of coefficients
        if (size == 0) return;
        if (size < PackSize) {
            visitor.bitpack(xpr.bitpack(0), size);
            return;
        }
        std::size_t k = 0, i = 0;   // k: current bitpack, i: maximum coefficient index processed
        for (; i + PackSize - 1 < size; i += PackSize) {
            visitor.bitpack(xpr.bitpack(k));
	    if(visitor) return;
            k++;
        }
        if (i < size) visitor.bitpack(xpr.bitpack(k), size - i);
        return;
    };
};

// evaluates true if all coefficients in the binary expression are true
template <typename XprType> struct bitpack_all_visitor {
    using BitPackType = typename XprType::BitPackType;
    bool res = true;
    inline void bitpack(BitPackType b) { res &= (~b == 0); }
    inline void bitpack(BitPackType b, std::size_t size) { res &= (~((((BitPackType)1 << size) - 1) | b) == 0); }
    operator bool() const { return res == false; }   // stop if already false
};

// evaluates true if at least one coefficient in the binary expression is true
template <typename XprType> struct bitpack_any_visitor {
    using BitPackType = typename XprType::BitPackType;
    bool res = false;
    inline void bitpack(BitPackType b) { res |= (b != 0); }
    inline void bitpack(BitPackType b, std::size_t size) { res |= (((~(BitPackType)0 >> size) & b) != 0); }
    operator bool() const { return res == true; }   // stop if already true
};

template <int Rows, int Cols, typename XprType> class BinMtxBase {
   protected:
    std::size_t n_rows_ = 0, n_cols_ = 0;
    std::size_t n_bitpacks_ = 0;   // number of required bitpacks
   public:
    using BitPackType = std::uintmax_t;
    static constexpr std::size_t PackSize = sizeof(BitPackType) * 8; // number of bits in a packet
  
    BinMtxBase() = default;
    BinMtxBase(std::size_t n_rows, std::size_t n_cols) :
        n_rows_(n_rows), n_cols_(n_cols), n_bitpacks_(1 + std::ceil((n_rows_ * n_cols_) / PackSize)) {};
    // getters
    inline std::size_t rows() const { return n_rows_; }
    inline std::size_t cols() const { return n_cols_; }
    inline std::size_t bitpacks() const { return n_bitpacks_; }
    inline std::size_t size() const { return n_rows_ * n_cols_; }
    XprType& get() { return static_cast<XprType&>(*this); }
    const XprType& get() const { return static_cast<const XprType&>(*this); }
    // access operator on base type E
    bool operator()(std::size_t i, std::size_t j) const {
        fdapde_assert(i < n_rows_ && j < n_cols_);
        return get().operator()(i, j);
    }
    // access to i-th bitpack of the expression
    BitPackType bitpack(std::size_t i) const { return get().bitpack(i); }
    // send matrix to out stream
    friend std::ostream& operator<<(std::ostream& out, const BinMtxBase& m) {
      // assign to temporary (triggers fast bitwise evaluation)
      BinaryMatrix<Rows, Cols> tmp;
      tmp = m;
      for (std::size_t i = 0; i < tmp.rows() - 1; ++i) {
            for (std::size_t j = 0; j < tmp.cols(); ++j) { out << tmp(i, j); }
            out << "\n";
        }
        // print last row without carriage return
        for (std::size_t j = 0; j < tmp.cols(); ++j) { out << tmp(tmp.rows() - 1, j); }
        return out;
    }
    // expression bitwise NOT
    BinMtxNegationOp<Rows, Cols, XprType> operator~() const {
        return BinMtxNegationOp<Rows, Cols, XprType>(get(), n_rows_, n_cols_);
    }
    // block-type indexing
    BinMtxBlock<1, Cols, XprType> row(std::size_t row) { return BinMtxBlock<1, Cols, XprType>(get(), row); }
    BinMtxBlock<1, Cols, const XprType> row(std::size_t row) const {
        return BinMtxBlock<1, Cols, const XprType>(get(), row);
    }
    BinMtxBlock<Rows, 1, XprType> col(std::size_t col) { return BinMtxBlock<Rows, 1, XprType>(get(), col); }
    BinMtxBlock<Rows, 1, const XprType> col(std::size_t col) const {
        return BinMtxBlock<Rows, 1, const XprType>(get(), col);
    }
    template <int Rows_, int Cols_>
    BinMtxBlock<Rows_, Cols_, XprType> block(std::size_t start_row, std::size_t start_col) {
        return BinMtxBlock<Rows_, Cols_, XprType>(get(), start_row, start_col);
    }
    BinMtxBlock<Dynamic, Dynamic, XprType>
    block(std::size_t start_row, std::size_t start_col, std::size_t block_rows, std::size_t block_cols) {
        return BinMtxBlock<Dynamic, Dynamic, XprType>(get(), start_row, start_col, block_rows, block_cols);
    }
    // visitors support
    inline bool all() const { return visit_apply_<bitpack_all_visitor<XprType>>(); }
    inline bool any() const { return visit_apply_<bitpack_any_visitor<XprType>>(); }
   private:
    template <typename Visitor> inline bool visit_apply_() const {
        Visitor visitor;
        bitpack_visit<XprType, Visitor>::run(get(), visitor);
        return visitor.res;
    }
};

// alias export for binary vectors
template <int Rows> using BinaryVector = BinaryMatrix<Rows, 1>;
  
}   // namespace core
}   // namespace fdapde


#endif   // __BINARY_MATRIX_H__
