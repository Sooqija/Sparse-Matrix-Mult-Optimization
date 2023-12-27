#include <iostream>
#include <vector>
#include <random>

class Complex {
 private:
    double rl;
    double im;
 public:
    explicit Complex(double _rl = 0, double _im = 0): rl(_rl), im(_im) {}
    Complex(const Complex& Tmp): rl(Tmp.rl), im(Tmp.im) {}
    Complex& operator=(Complex Tmp) {
        this->rl = Tmp.rl;
        this->im = Tmp.im;
        return *this;
    }
    std::vector<Complex> InitVec(std::vector<double> rls = std::vector<double>(),
                                 std::vector<double> ims = std::vector<double>()) {
        std::vector<Complex> Ent(rls.size());
        if (rls.size() == ims.size()) {
            for (size_t i = 0; i < rls.size(); i++) {
                Ent[i] = Complex(rls[i], ims[i]);
            }
        } else {
            for (size_t i = 0; i < rls.size(); i++) {
                Ent[i] = Complex(rls[i], 0);
            }
        }
        return Ent;
    }
    double GetRl() { return this->rl; }
    double GetIm() { return this->im; }
    void SetRl(double tmp) { this->rl = tmp; }
    void SetIm(double tmp) { this->im = tmp; }
    Complex operator+(Complex Tmp) {
        Complex Ans;
        Ans.rl = this->rl + Tmp.rl;
        Ans.im = this->im + Tmp.im;
        return Ans;
    }
    Complex& operator+=(Complex Tmp) {
        this->rl += Tmp.rl;
        this->im += Tmp.im;
        return *this;
    }
    Complex operator*(Complex Tmp) {
        Complex Ans;
        Ans.rl = this->rl * Tmp.rl - this->im * Tmp.im;
        Ans.im = this->rl * Tmp.im + this->im * Tmp.rl;
        return Ans;
    }
    bool operator==(Complex Tmp) {
        if ((this->rl - Tmp.rl < 0.00001) && (this->im - Tmp.im < 0.00001)) {
            return true;
        } else {
            return false;
        }
    }
    bool operator!=(Complex Tmp) { return !(*this == Tmp); }
    bool IsNotZero() {
        bool ans = false;
        const double ZeroLike = 0.000001;
        if ((fabs(this->rl) > ZeroLike) || (fabs(this->im) > ZeroLike)) {
            ans = true;
        }
        return ans;
    }
    ~Complex() {}
};

class Matrix {
 private:
    int size;
    int nonzv;
    std::vector<Complex> values;
    std::vector<int> irows;
    std::vector<int> shiftcols;

 public:
    Matrix(int _size = 0, int _non = 0, const std::vector<Complex>& _Entry =
           std::vector<Complex>(), const std::vector<int>& _irows =
           std::vector<int>(), const std::vector<int>& _shtcols =
           std::vector<int>()): size(_size), nonzv(_non), values(_Entry),
           irows(_irows), shiftcols(_shtcols) {}
    Matrix(const Matrix& Tmp): size(Tmp.size), nonzv(Tmp.nonzv), values(Tmp.values),
                               irows(Tmp.irows), shiftcols(Tmp.shiftcols) {}
    Matrix& operator=(const Matrix& Tmp) {
        this->size = Tmp.size;
        this->nonzv = Tmp.nonzv;
        this->values = Tmp.values;
        this->irows = Tmp.irows;
        this->shiftcols = Tmp.shiftcols;

        return *this;
    }
    int GetSize() { return this->size; }
    int GetNon() { return this-> nonzv; }
    std::vector<Complex> GetEntry() { return this->values; }
    std::vector<int> GetIrows() { return this->irows; }
    std::vector<int> GetShtcols() { return this->shiftcols; }
    Matrix& ClearMatrix() {
        this->size = 0;
        this->nonzv = 0;
        this->values.clear();
        this->irows.clear();
        this->shiftcols.clear();

        return *this;
    }
    Matrix& RandomMatrix(int size, int dist, int cnt = -1, int seed = 0) {
        this->ClearMatrix();
        this->size = size;
        std::mt19937 gen(time(0));
        gen.seed(seed);
        if (cnt < 0) {
            cnt = static_cast<int>(size * 0.01);
        }
        this->nonzv = cnt * size;
        this->values.resize(cnt * size);
        this->irows.resize(cnt * size);
        this->shiftcols.resize(size + 1);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < cnt; j++) {
                bool flag;
                do {
                    this->irows[i * cnt + j] = gen() % size + 1;
                    flag = false;
                    for (int k = 0; k < j; k++) {
                        if (this->irows[i * cnt + j] == this->irows[i * cnt + k]) {
                            flag = true;
                        }
                    }
                } while (flag);
            }
            for (int j = 0; j < cnt - 1; j++) {
                for (int k = 0; k < cnt - 1; k++) {
                    if (this->irows[i * cnt + k] > this->irows[i * cnt + k]) {
                        int tmp = this->irows[i * cnt + k];
                        this->irows[i * cnt + k] = this->irows[i * cnt + k + 1];
                        this->irows[i * cnt + k] = tmp;
                    }
                }
            }
        }
        for (int i = 0; i < cnt * size; i++) {
            this->values[i].SetRl(gen() % dist + 1);
            this->values[i].SetIm(gen() % dist + 1);
        }
        int sum = 1;
        for (int i = 0; i < size + 1; i++) {
            this->shiftcols[i] = sum;
            sum += cnt;
        }

        return *this;
    }
    bool operator==(Matrix Tmp) {
            bool ans = true;
        if (this->nonzv != Tmp.nonzv) {
            return false;
        }
        if (this->size != Tmp.size) {
            return false;
        }
        for (int i = 0; i < this->nonzv; i++) {
            if (this->values[i] != Tmp.values[i]) {
                return false;
            }
        }
        if (this->irows != Tmp.irows) {
            return false;
        }
        if (this->shiftcols != Tmp.shiftcols) {
            return false;
        }
        return ans;
    }
    bool operator!=(const Matrix& Tmp) { return !(*this == Tmp); }
    Matrix T() {
        Matrix Ans;
        Ans.nonzv = this->nonzv;
        Ans.size = this->size;
        Ans.values.resize(this->nonzv);
        Ans.irows.resize(this->nonzv);
        Ans.shiftcols.resize(this->size + 1);

        for (int i = 0; i < this->nonzv; i++) {
            Ans.shiftcols[this->irows[i] - 1]++;
        }
        int sum = 1;
        for (int i = 0; i < this->size + 1; i++) {
            int tmp = Ans.shiftcols[i];
            Ans.shiftcols[i] = sum;
            sum += tmp;
        }
        std::vector<int> shiftcols_tmp = Ans.shiftcols;
        for (int i = 0; i < this->size; i++) {
            for (int j = this->shiftcols[i]; j < this->shiftcols[i + 1]; j++) {
                int r_index = this->irows[j - 1];
                int i_index = shiftcols_tmp[r_index - 1];
                Ans.irows[i_index - 1] = i + 1;
                Ans.values[i_index - 1] = this->values[j - 1];
                shiftcols_tmp[r_index - 1]++;
            }
        }
        return Ans;
    }
    Matrix operator*(Matrix B) {
        Matrix A = this->T();

        std::vector<Complex> values;
        std::vector<int> irows_res;
        std::vector<int> shiftcols_res = { 1 };

        int nonzv_counter = 1;
        for (int j = 0; j < B.size; j++) {
            std::vector<int> ip(A.size, 0);
            for (int i = B.shiftcols[j]; i < B.shiftcols[j+1]; i++) {
                ip[B.irows[i-1]-1] = i;
            }
            for (int i = 0; i < A.size; i++) {
                Complex sum;
                for (int k = A.shiftcols[i]; k < A.shiftcols[i+1]; k++) {
                    int p = ip[A.irows[k-1]-1];
                    if (p) {
                        sum += B.values[p-1] * A.values[k-1];
                    }
                }
                if (sum.IsNotZero()) {
                        irows_res.push_back(i+1);
                        values.push_back(sum);
                        nonzv_counter++;
                }
            }
            shiftcols_res.push_back(nonzv_counter);
        }
        Matrix Ans(A.size, nonzv_counter-1, values, irows_res, shiftcols_res);
        return Ans;
    }
    ~Matrix() {
        this->values.clear();
        this->irows.clear();
        this->shiftcols.clear();
    }
};



int main()
{
    int size = 10000;
    int dist = 1000;
    int cnt = 50;
    Matrix A;
    A.RandomMatrix(size, dist, cnt, 0);
    Matrix B;
    B.RandomMatrix(size, dist, cnt, 1);

    printf("Mult started\n");
    A * B;
    printf("Mult completed\n");
}
