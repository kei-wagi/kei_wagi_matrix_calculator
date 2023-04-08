/*
 * 「結果」が行列である電卓. 
 * 行列クラスを定義し, とりあえず演算としては加算や単位行列等を定義している. 
 * コンパイル & 実行：
 * javac Calculator.java IntCalc.java MemoCalc.java MatrixCalc.java
 * java MatrixCalc
 */

import java.util.*;


import java.io.*;
import java.math.*;

/**
 * 電卓の「結果」として使う行列を表すクラス. 
 * 行列の要素は {@code double} の2次元配列で保持する. 
 * 加算や単位行列生成などの演算や, 行列を文字列から読み込む機能を提供する. 
 */
class Matrix {
    /**
     * 行列の行数. 
     */
    final int m;
    /**
     * 行列の列数. 
     */
    final int n;
    /**
     * 行列の要素. 
     * 並びは自然な並びで： {@code vals[i][j]} が (i, j) 要素. 
     */
    double [][] vals;
    /**
     * {@code m}×{@code n} のゼロ行列を作るコンストラクタ. 
     * @param m 行数 
     * @param n 列数 
     */
    Matrix(int m, int n) {
        this.m = m;
        this.n = n;
        vals = new double[m][n];
    }
    /**
     * 与えられた行列をコピーするコンストラクタ. 
     * @param mat コピー元の行列. 
     */
    Matrix(Matrix mat) {
        this(mat.m, mat.n);
        copy(mat.vals);
    }
    /**
     * 与えられた2次元配列の内容を自身の要素としてコピーする. 
     * 次元は矛盾しないとする（与えられた2次元配列の方が大きければ良い）. 
     * @param vals コピー元の2次元配列. 
     */
    void copy(double [][] vals) {
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                this.vals[i][j] = vals[i][j];
            }
        }
    }

    public static Matrix buffer_LU;
    public static Matrix buffer_Eigen;

    /**
     * 与えら得た行列がサイズ違いで自身に加減算できないときに {@code true} を返す. 
     * @param mat 行列. 
     * @return {@code mat} と自身のサイズが同じでないときに {@code true}. 
     */
    boolean sizeMismatch(Matrix mat) {
        return (mat.m != m) || (mat.n != n);
    }
    
    
    
    
    /**
     * 与えら得た行列がサイズ違いで自身に乗算できないときに {@code true} を返す. 
     * @param mat 行列. 
     * @return {@code mat} と自身のサイズが同じでないときに {@code true}. 
     */
    boolean sizeMismatch_mul(Matrix mat) {
        return (mat.m != this.n && this.n !=1 && this.m !=1 );
    }



    /**
     * 与えられた行列と自身の加算結果の行列を新たに生成して返す. 
     * @param mat 加算する行列
     * @return 行列加算 {@code this} + {@code mat} の結果となる行列. 
     *         サイズ違いなどで計算不可能な場合には {@code null}. 
     */

    Matrix add(Matrix mat) {
        // 計算できないときには null を返す. 
        if(mat == null || sizeMismatch(mat)) return null;
        // あとは単純な加算
        Matrix ret = new Matrix(m, n);
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                ret.vals[i][j] = this.vals[i][j] + mat.vals[i][j];
            }
        }
        return ret;
    }

    /**
     * 与えられた行列と自身の減算結果の行列を新たに生成して返す. 
     * @param mat 減算する行列
     * @return 行列減算 {@code this} - {@code mat} の結果となる行列. 
     *         サイズ違いなどで計算不可能な場合には {@code null}. 
     */
    Matrix sub(Matrix mat) {
        // 計算できないときには null を返す. 
        if(mat == null || sizeMismatch(mat)) return null;
        // あとは単純な減算
        Matrix ret = new Matrix(m, n);
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                ret.vals[i][j] = this.vals[i][j] - mat.vals[i][j];
            }
        }
        return ret;
    }


    Matrix mul(Matrix mat) {
        // 計算できないときには null を返す. 
        if(mat == null || sizeMismatch_mul(mat)) return null;
        // あとは単純な乗算
        Matrix ret = new Matrix(this.m, mat.n);
        Matrix ret1 = new Matrix(mat.m, mat.n);
        if(this.n ==1 && this.m == 1){
            for(int i = 0; i < mat.m; i++) {
                for(int j = 0; j < mat.n; j++) {
                    ret1.vals[i][j] = this.vals[0][0] * mat.vals[i][j];
                }
            }
            ret = ret1;
        }else{
            for(int i = 0; i < this.m; i++) {
                for(int j = 0; j < mat.n; j++) {
                    for(int k = 0; k < this.n; k++){
                        ret.vals[i][j] += this.vals[i][k] * mat.vals[k][j];
                    }
                }
            }
        }

        return ret;
    }



    public static Matrix LU(Matrix mat, boolean a)
    {
        if(mat == null || Matrix.Eigen(mat, true).vals[0][0] == 0 || mat.m != mat.n || mat.vals[0][0] == 0 ) return null;
        Matrix L = new Matrix(mat.m, mat.m);
        Matrix U = new Matrix(mat.m, mat.m);
        int i, j, k; 
        double T;
       L.vals[0][0] = 1;
       for(j = 0; j < mat.m; j++){
            U.vals[0][j] = mat.vals[0][j];
       }
       for(i = 1; i < mat.m; i++){
            U.vals[i][0] = L.vals[0][i] = 0;
            L.vals[i][0] = mat.vals[i][0] / U.vals[0][0];
       }
       for(i = 1; i < mat.m; i++){
            L.vals[i][i] = 1; T = mat.vals[i][i];
            for(k=0; k<=i; k++){
            T -= L.vals[i][k] * U.vals[k][i];
            }
            U.vals[i][i] = T;
            for(j=i+1; j < mat.m; j++){
                U.vals[j][i] = L.vals[i][j] = 0;
                T = mat.vals[j][i];
                for(k = 0; k <= i; k++){
                    T -= L.vals[j][k] * U.vals[k][i];
                }
                L.vals[j][i] = T / U.vals[i][i];
                T = mat.vals[i][j];
                for(k = 0; k <= i; k++){
                    T -= L.vals[i][k] * U.vals[k][j];  
                } 
                U.vals[i][j] = T;
            }
        }
        if(a){
            buffer_LU = U;
            return L;
        }else{
            buffer_LU = L;
            return U;
        }
    }


    public static Matrix Inv(Matrix mat)
    {   if(Matrix.det(mat).vals[0][0] == 0){
            return null;
        }
        int i,j,k,l;
        int index;
        double Alfa, tmp;
        Matrix mat_1 = new Matrix(mat.m, mat.m);
        int N2 = mat.m * 2;    
/*         for(i=0;i<mat.m;i++){
            for(j=0;j<mat.m;j++){
                mat_1.vals[i][j] = (i==j ? 1 : 0);  
            } 
        } */
        mat_1 = eye(mat.m);

        for(k = 0; k < mat.m; k++){
            if(mat.vals[k][k] == 0){
                index = k;
                double max = mat.vals[k][k];
                for(l=0;l<mat.m;l++){
                    if(max < mat.vals[l][k]){
                        max = mat.vals[l][k];
                        index = l;
                    }
                }
                if(max == 0){
                    return null;
                }else{
                    for(l=0;l<mat.m;l++){
                        tmp = mat.vals[index][l];
                        mat.vals[index][l] = mat.vals[k][l];
                        mat.vals[k][l] = tmp;
                        tmp = mat_1.vals[index][l];
                        mat_1.vals[index][l] = mat_1.vals[k][l];
                        mat_1.vals[k][l] = tmp;
                    }
                }
                Alfa=1 / mat.vals[k][k];

            }else{
                Alfa=1 / mat.vals[k][k];
            }       
            int k1=k+1;
            for(j=k1;j<N2;j++){
                if(j < mat.m){
                   mat.vals[k][j] *= Alfa; 
                }else{
                    mat_1.vals[k][j - mat.m] *= Alfa;
                }
                
            }
            for(i=0;i<mat.m;i++){
                if(i != k){
                    Alfa=mat.vals[i][k];
                    for(j=k1;j<N2;j++){
                        if(j < mat.m){
                           mat.vals[i][j] -= Alfa* mat.vals[k][j];  
                        }else{
                            mat_1.vals[i][j-mat.m] -= Alfa* mat_1.vals[k][j-mat.m]; 
                        }
                         
                    }
                }
            }
        } 
        return mat_1;
    }

	private static interface Access {
		
		int size();
		double get(int r, int c);
	}
	
	private static class WrapEntity implements Access {

        private Access ma;
		private int row;
		private int col;
		private int size;
		
		public WrapEntity(Access mat, int r, int c){
            this.ma = mat;
			this.row = r;
			this.col = c;
			this.size = mat.size() - 1;
		}
		
		public double get(int r, int c) {
			r = r>=row ? r+1: r;
			c = c>=col ? c+1: c;
			return ma.get(r, c);
		}

		public int size() {
			return size;
		}
		
	}
	
	private static class Entity implements Access {
		
		private Matrix ma;
		private int row;		
		private int col;
		private int size;
		
		public Entity(Matrix mat, int r, int c){
            this.ma = mat;
			this.row = r;
			this.col = c;
			this.size = mat.m - 1;
		}
		
		public int size(){
			return size;
		}
		
		public double get(int r, int c){
			r = r>=row ? r+1: r;
			c = c>=col ? c+1: c;
			return ma.vals[r][c];
		}
	}
	
	private static double sub1(Access mat){
		if(mat.size() == 1){
			return mat.get(0, 0);
		}
		if(mat.size()==2){
			return mat.get(0, 0)* mat.get(1, 1)-mat.get(0, 1)*mat.get(1, 0);
		}
		// 行列式の計算
		double det = 0;
		int length = mat.size();
		for(int i=0;i<length;i++){
			double v = mat.get(0, i);
			if(v!=0){
				double d = v*sub1(new WrapEntity(mat,0 ,i));
				det += i%2==0? d: -d;
			}
		}
		return det;
	}
	
	public static Matrix det(Matrix mat){
		// 引数チェック
		if(mat.n!=mat.m){
			return null;
        }

		// 行列式の計算
		double det = 0;
		for(int i=0;i<mat.m;i++){
			double v = mat.vals[0][i];
			if(v!=0){
				double d = v*sub1(new Entity(mat,0,i));
				det += i%2==0? d: -d;
			}
		}
        Matrix answer = eye(1);
        answer.vals[0][0]=det;
		return answer;
    }



    public static Matrix Eigen(Matrix mat, boolean a){
        double lambda = 0;
        double EPS = 0.0000001;
        double T;
        int i,j, K1, K2, KK, IterMax;
        Matrix Mtemp = new Matrix(2, mat.m);
        K1 = 1;
        K2 = 0;
        IterMax = 100;
        Mtemp.vals[K2][0] = 1; 
        
        double E = EPS * 100;
        int Iter = 0;
        while(Iter < IterMax && E > EPS){
            Iter++; 
            KK = K2; 
            K2 = K1; 
            K1 = KK; 
            for(i = 0; i < mat.m; i++){
                T = 0;
                for(j = 0; j < mat.m; j++){
                    T += mat.vals[i][j] * Mtemp.vals[K1][j]; 
                } 
                Mtemp.vals[K2][i] = T;
            }
 
            T = 0; 
            for(i = 0; i < mat.m; i++){
                T += Mtemp.vals[K2][i]*Mtemp.vals[K2][i];
            } 
            lambda = Math.sqrt(T);
            if (Math.abs(lambda ) < EPS){
                E = EPS * 100; 
                break;
            }
            for(i = 0; i < mat.m; i++){
                Mtemp.vals[K2][i] = Mtemp.vals[K2][i] / lambda;
            }
            double T1 = 0;
            double T2 = 0;
            for(i = 0; i < mat.m; i++){  
                double DX = Mtemp.vals[K2][i] - Mtemp.vals[K1][i];
                T1 += DX * DX; 
                T2 += Mtemp.vals[K2][i] * Mtemp.vals[K2][i];
            }
            if(Math.abs(T2) < EPS){  
                E = EPS * 100; 
                break;
            }
            E = T1 / T2;
       }
       Matrix value = eye(1);
       Matrix vector = new Matrix(mat.m, 1); ;
       value.vals[0][0]=lambda;
       for(i = 0; i < mat.m; i++){
            vector.vals[i][0] = Mtemp.vals[K2][i];
       }
       if(a){
            buffer_Eigen = vector;
            return value;
        }else{
            buffer_Eigen = value;
            return vector;
        }


    }



    Matrix Sys_Eq(Matrix mat){
        if(mat == null || this.m != this.n || this.m != mat.m) return null;
        int i,k; 
        double T;
        Matrix L = Matrix.LU(this,true);
        Matrix U = buffer_LU;
        int N = mat.m;

        Matrix Y = new Matrix(N, 1);
        Matrix X = new Matrix(N, 1);
        
        for(i=0;i<N;i++){ 
            T = mat.vals[i][0]; 
            for(k=0;k<=i-1;k++){
                T -= L.vals[i][k]*Y.vals[k][0];
            }
            Y.vals[i][0]=T;
        }
        for(i=N-1;i>=0;i--){
            T = Y.vals[i][0];
            for(k=i+1;k<N;k++){
                T -= U.vals[i][k]*X.vals[k][0];
            }
            X.vals[i][0]=T/U.vals[i][i];          
        }
        return X;
    }




    /**
     * 与えられたサイズの単位行列を新たに生成して返す. 
     * @param n 生成する行列のサイズ
     * @return {@code n}×{@code} の単位行列
     */
    public static Matrix eye(int n) {
        Matrix ret = new Matrix(n, n); // これは nxn のゼロ行列
        for(int i = 0; i < n; i++) {
            ret.vals[i][i] = 1;  // 対角に 1 を入れる
        }
        return ret;
    }
    /**
     * 「ブロック」から与えられたサイズの単位行列を新たに生成して返す. 
     * @param block 電卓から受け取る「ブロック」.
     * @return {@code n}×{@code} の単位行列. 行のサイズの食い違いなどで生成に失敗したら {@code null}
     */
    public static Matrix read(final List<String> block) {
        try {
            int m = block.size() - 1;  // 一行目は行列の中身ではないので無視して行数を決める
            int n = -1;                // 列数は, 最初の行を見て決める
            Matrix ret = null;
            for(int i = 0; i < m; i++) {
                StringTokenizer st = new StringTokenizer(block.get(i+1));
                ArrayList<String> vs = new ArrayList<String>();
                while(st.hasMoreTokens()) vs.add(st.nextToken());
                if(n < 0) {  // 最初の行でサイズ確定 → 決定したサイズの行列をここで生成
                    n = vs.size();
                    ret = new Matrix(m, n);
                } else if(n != vs.size()) {// 行の間でサイズの食い違いがあったら null
                    return null;
                }
                // 要素をコピー
                int j = 0;
                for(String s : vs) {
                    ret.vals[i][j++] = Double.parseDouble(s);
                }
            }
            return ret;
        } catch(Exception e) { // なにか変な例外が生じた際にも生成失敗
        }
        return null;
    }
    /**
     * 行列を表す文字列を返す. 
     * 例えば次のような文字列となる. 
     * <p><blockquote><pre>{@code
     * [   2.000    3.000    4.000]
     * [   5.000    6.000    7.000]
     * }</pre></blockquote><p>
     * （各行の開始と終わりに [ と ] が置かれ, 各要素は固定幅（8.3f）で表示. 
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        for(int i = 0; i < m; i++) {
            sb.append("[");
            for(int j = 0; j < n; j++) {
                if(j > 0) sb.append(" ");
                sb.append(String.format("%1$8.3f", vals[i][j]));
            }
            sb.append("]");
            if(i < m - 1) sb.append("\n");
        }
        return sb.toString();
    }
}

/**
 * 行列加算を入力して現在の「結果」をその行列にする「コマンド」. 
 * <p><blockquote><pre>{@code
 * mat :
 *  TAB  a_11 ... a_1m
 *  TAB  a_21 ... a_2m
 *    ...
 *  TAB  a_n1 ... a_nm
 * }</pre></blockquote><p>
 * という, 1行目が {@code mat} である複数行「ブロック」を受け付け, 入力された行列を「結果」として返す. 
 * 行列の入力は, 各行の要素を TAB 始まりの各行に空白区切りで並べて入力する. 
 * 例えば, 次のような「ブロック」を入力として受け付ける（2行目以降は TAB を先頭に入力すること）. 
 * <p><blockquote><pre>{@code
 * mat :
 *      2 3 4
 *      5 6 7
 * }</pre></blockquote><p>
 */
class MatrixValue implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix r) {
        if(block.size() <= 1) return null;
        if(ts.length == 1 && "mat".equals(ts[0])) {
            // 実際の読み込みは Matrix クラスに任せる
            return Matrix.read(block);
        }
        return null;
    }
}

/**
 * 単位行列を現在の「結果」にする「コマンド」. 
 * {@code eye} の後に整数値が並ぶ 1行の「ブロック」を受け付け, その整数値のサイズの単位行列を「結果」として返す. 
 * 例えば, 次のような「ブロック」を入力として受け付ける（2x2 の単位行列になる）. 
 * <p><blockquote><pre>{@code
 * eye 2
 * }</pre></blockquote><p>
 */
class IdentityMatrix implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix r) {
        if(block.size() != 1) return null;
        if(ts.length == 2 && "eye".equals(ts[0])) {
            // 単位行列の実際の生成は Matrix クラスにまかせる
            return Matrix.eye(Integer.parseInt(ts[1]));
        }
        return null;
    }
}

/**
 * ゼロ行列を現在の「結果」にする「コマンド」. 
 * <p><blockquote><pre>{@code
 * zero n
 * }</pre></blockquote><p>
 * のような 1行「ブロック」を受け付け, その整数値 {@code n} のサイズのゼロ行列を「結果」として返す. 
 * 例えば, 次のような「ブロック」を入力として受け付ける（2x2 のゼロ行列になる）. 
 * <p><blockquote><pre>{@code
 * zero 2
 * }</pre></blockquote><p>
 */
class ZeroMatrix implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix r) {
        if(block.size() != 1) return null;
        if(ts.length == 2 && "zero".equals(ts[0])) {
            int n = Integer.parseInt(ts[1]);
            return new Matrix(n, n); // コンストラクタで生成した状態がゼロ行列なので
        }
        return null;
    }
}

/**
 * 行列加算の「コマンド」. 
 * <p><blockquote><pre>{@code
 * add :
 *  TAB  a_11 ... a_1m
 *  TAB  a_21 ... a_2m
 *    ...
 *  TAB  a_n1 ... a_nm
 * }</pre></blockquote><p>
 * という, 1行目が {@code add} である複数行「ブロック」を受け付け, 入力された行列を足した「結果」を返す. 
 * 行列の入力は, 各行の要素を TAB 始まりの各行に空白区切りで並べて入力する. 
 * 例えば, 次のような「ブロック」を入力として受け付ける（2行目以降は TAB を先頭に入力すること）. 
 * <p><blockquote><pre>{@code
 * add :
 *      1 0 1
 *      0 1 0
 * }</pre></blockquote><p>
 * <br />
 * もしくは, {@code add} の後ろに変数名を書いた 1行の「ブロック」を受け付け, 
 * 変数に保存された行列を足した「結果」を返す.  
 */
class MatrixAdd extends CommandWithMemory<Matrix> {
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    MatrixAdd(Memory<Matrix> mem) {
        super(mem); // 親のコンストラクタをそのまま呼ぶだけ
    }
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        // 行列の値を直接書く場合
        if(block.size() > 1 && ts.length == 1 && "add".equals(ts[0])){
            // 実際の読み込みと加算は Matrix クラスに任せる
            Matrix v = Matrix.read(block);
            return res.add(v);
        }
        // 行列を保存した変数が指定された場合
        if(block.size() == 1 && ts.length == 2 && "add".equals(ts[0])) {
            // 変数の値をメモリから取得
            Matrix v = mem.get(ts[1]);
            return res.add(v); // 実際の加算は Matrix クラス任せ
        }
        return null;
    }
}











class MatrixSub extends CommandWithMemory<Matrix> {
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    MatrixSub(Memory<Matrix> mem) {
        super(mem); // 親のコンストラクタをそのまま呼ぶだけ
    }
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        // 行列の値を直接書く場合
        if(block.size() > 1 && ts.length == 1 && "sub".equals(ts[0])){
            // 実際の読み込みと減算は Matrix クラスに任せる
            Matrix v = Matrix.read(block);
            return res.sub(v);
        }
        // 行列を保存した変数が指定された場合
        if(block.size() == 1 && ts.length == 2 && "sub".equals(ts[0])) {
            // 変数の値をメモリから取得
            Matrix v = mem.get(ts[1]);
            return res.sub(v); // 実際の減算は Matrix クラス任せ
        }
        return null;
    }
}

class MatrixMul extends CommandWithMemory<Matrix> {
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    MatrixMul(Memory<Matrix> mem) {
        super(mem); // 親のコンストラクタをそのまま呼ぶだけ
    }
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        // 行列の値を直接書く場合
        if(block.size() > 1 && ts.length == 1 && "mul".equals(ts[0])){
            // 実際の読み込みと乗算は Matrix クラスに任せる
            Matrix v = Matrix.read(block);
            return res.mul(v);
        }
        // 行列を保存した変数が指定された場合
        if(block.size() == 1 && ts.length == 2 && "mul".equals(ts[0])) {
            // 変数の値をメモリから取得
            Matrix v = mem.get(ts[1]);
            return res.mul(v); // 実際の乗算は Matrix クラス任せ
        }
        return null;
    }
}



class MatrixLU implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        if(block.size() > 1) return null;
        if(ts.length == 1 && "lu_l".equals(ts[0])){
            return Matrix.LU(res, true);
        }else if(ts.length == 1 && "lu_u".equals(ts[0])){
            return Matrix.LU(res, false);
        }
        return null;
    }
}


class MatrixInv implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        if(block.size() > 1) return null;
        if(ts.length == 1 && "inv".equals(ts[0])){
            return Matrix.Inv(res);
        }
        return null;
    }
}

class MatrixDet implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        if(block.size() > 1) return null;
        if(ts.length == 1 && "det".equals(ts[0])){
            return Matrix.det(res);
        }
        return null;
    }
}

class MatrixEigen implements Command<Matrix> {
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        if(block.size() > 1) return null;
        if(ts.length == 1 && "eigen_val".equals(ts[0])){
            return Matrix.Eigen(res, true);
        }else if(ts.length == 1 && "eigen_vec".equals(ts[0])){
            return Matrix.Eigen(res, false);
        }
        return null;
    }
}




class MatrixSys_eq extends CommandWithMemory<Matrix> {
    MatrixSys_eq(Memory<Matrix> mem) {
        super(mem);
    }
    public Matrix tryExec(final String [] ts, final List<String> block, final Matrix res) {
        if(block.size() > 1 && ts.length == 1 && "sys_eq".equals(ts[0])){
            Matrix v = Matrix.read(block);
            return res.Sys_Eq(v);
        }
        if(block.size() == 1 && ts.length == 2 && "sys_eq".equals(ts[0])) {
            Matrix v = mem.get(ts[1]);
            return res.Sys_Eq(v); 
        }
        return null;
    }
}

/**
 * 行列電卓を作成して動作させるクラス. 
 * 例えば, ターミナルで次のような実行ができる. 
 * <p><blockquote><pre>{@code
 * $ javac Calculator.java IntCalc.java MemoCalc.java MatrixCalc.java
 * $ java MatrixCalc
 * [   0.000    0.000]
 * [   0.000    0.000]
 * >> mat :
 * ..      2 3 4
 * ..      5 6 7
 * ..
 * [   2.000    3.000    4.000]
 * [   5.000    6.000    7.000]
 * >> store x
 * [   2.000    3.000    4.000]
 * [   5.000    6.000    7.000]
 * >> add :
 * ..      1 0 1
 * ..      0 1 0
 * ..
 * [   3.000    3.000    5.000]
 * [   5.000    7.000    7.000]
 * >> add x
 * [   5.000    6.000    9.000]
 * [  10.000   13.000   14.000]
 * >> eye 2
 * [   1.000    0.000]
 * [   0.000    1.000]
 * >>
 * }</pre></blockquote><p>
 * これは, 最初に「結果」が 2x2 のゼロ行列で電卓が動き始め, 
 * まずは最初のプロンプトの後ろで {@code mat : } と打って複数行の「ブロック」の入力を始め, 
 * 続く TAB 始まりの 2行に 2x3 行列の各行の要素を書き, 続く空行で「ブロック」の入力を終え, 
 * その「ブロック」に対して {@code MatrixValue} が実行されて「結果」がその行列になり, 
 * 続いて {@code store x} と入力し, 
 * {@code LoadStore} が動作してその行列が変数 x に保存され, 
 * 続いて {@code add :} からの数行で同様に 2x3 行列を入力し, 
 * それに対して {@code MatrixAdd} が動作して「結果」が行列の和になり, 
 * さらに {@cdoe add x} と入力し, 
 * それに対して {@code MatrixAdd} が動作して変数 x に保存した行列が加算され, 
 * 最後に {@code eye 2} と入力し, 
 * それに対して {@code IdentityMatrix} が動いて「結果」が 2x2 の単位行列となった. 
 */
class MatrixCalc {
    /**
     * 電卓を作って実行する. 
     */
    public static void main(String [] args) throws Exception {
        // 行列を記憶する変数のための Memory インスタンス
        Memory<Matrix> mem = new Memory<Matrix>();
        // コマンドリストの作成
        ArrayList<Command<Matrix>> comms = new ArrayList<Command<Matrix>>();
        comms.add(new EmptyCommand<Matrix>());
        comms.add(new MatrixValue());
        comms.add(new IdentityMatrix());
        comms.add(new ZeroMatrix());
        comms.add(new MatrixAdd(mem));
        comms.add(new MatrixSub(mem));
        comms.add(new MatrixMul(mem));
        comms.add(new MatrixLU());
        comms.add(new MatrixInv());
        comms.add(new MatrixDet());
        comms.add(new MatrixEigen());
        comms.add(new MatrixSys_eq(mem));
        comms.add(new LoadStore<Matrix>(mem));
        comms.add(mem);
        // 入力は標準入力から
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        // 電卓の生成と実行
        Calculator<Matrix> c = new Calculator<Matrix>(br, comms);
        // 初期値は 2x2 のゼロ行列
        c.run(new Matrix(2,2));
    }
}
