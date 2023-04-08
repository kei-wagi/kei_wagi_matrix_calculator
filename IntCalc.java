

import java.util.*;
import java.io.*;
import java.math.*;

/**
 * 新たな整数値を入力するための「コマンド」. 
 * 10進数の値のみを書いた「ブロック」を受け付けて, その値を「結果」として返す. 
 */
class IntValue implements Command<BigInteger> {
    /**
     * 「ブロック」が 1行のみで, かつ, その行が 10進数ひとつのみからなる場合に, その10進数の値を「結果」として返す. 
     */
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length != 1) return null;
        try {
            return new BigInteger(ts[0]); // とりあえず BigInteger に変換させる. 失敗したら例外が飛ぶ. 
        } catch(Exception e) {  // 変換失敗時の例外を受け取り, null を返す（＝ この「コマンド」は実行できない）
            return null;
        }
    }
}

/**
 * 加算のための「コマンド」. 
 * <p><blockquote><pre>{@code
 * + n
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, その値 {@code n} を現在の「結果」に足した「結果」を返す. 
 */
class IntAdd implements Command<BigInteger> {
    /**
     * 「ブロック」が 1行のみで, かつ, その行が "+" と10進数ひとつのみからなる場合に, 
     * その10進数の値を与えられた「結果」に足した「結果」を返す. 
     * 面倒なので 10進数以外が来ないことを期待した手抜き実装. 真面目にやるなら例外をキャッチして null を返す. 
     */
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length == 2 && "+".equals(ts[0])) {
            BigInteger v = new BigInteger(ts[1]);
            return res.add(v);
        } else {
            return null;
        }
    }
}

/**
 * 減算のための「コマンド」. 
 * <p><blockquote><pre>{@code
 * - n
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, その値 {@code n} を現在の「結果」から引いた「結果」を返す. 
 */
class IntSub implements Command<BigInteger> {
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length == 2 && "-".equals(ts[0])) {
            BigInteger v = new BigInteger(ts[1]);
            return res.subtract(v);
        } else {
            return null;
        }
    }
}

/**
 * 乗除算のための「コマンド」. 
 * <p><blockquote><pre>{@code
 * * n
 * }</pre></blockquote><p>
 * ないし
 * <p><blockquote><pre>{@code
 * / n
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, 
 * その値{@code n} を現在の「結果」に掛けた「結果」, 
 * ないし, その値{@code n} で現在の「結果」を割った「結果」をそれぞれ返す. 
 * ひとつオブジェクトで複数の演算に対応するための例. 
 */
class IntMulDiv implements Command<BigInteger> {
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length != 2) return null;

        // ひとつめのトークンが"*" なら掛け算をして返す
        if("*".equals(ts[0])) return res.multiply(new BigInteger(ts[1])); 

        // ひとつめのトークンが"/" なら割り算をして返す
        if("/".equals(ts[0])) return res.divide(new BigInteger(ts[1]));

        // それ以外は実行できないので null を返す
        return null;
    }
}

/**
 * 符号反転のための「コマンド」. 
 * <p><blockquote><pre>{@code
 * neg
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, 現在の「結果」を符号反転した「結果」を返す. 
 * 単項演算（例えば他に sqrt とか）の例. 
 */
class IntNeg implements Command<BigInteger> {
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length != 1) return null;
        if("neg".equals(ts[0])) return res.negate();
        return null;
    }
}

/**
 * 四則演算と符号反転をもった整数電卓を作成し動作させるクラス. 
 * {@code main} メソッドで「コマンド」のリストを作り, 
 * それと標準入力から読み込む {@code BufferedReader} とで
 * 新たな電卓を生成し, それを実行している. <br />
 * 例えば, ターミナルで次のような実行ができる. 
 * <p><blockquote><pre>{@code
 * $ javac Calculator.java IntCalc.java
 * $ java IntCalc
 * 0
 * >> 123
 * 123
 * >> + 45
 * 168
 * >> 
 * }</pre></blockquote><p>
 * これは, 最初に「結果」が 0 で電卓が動き始め, まずは最初のプロンプト（{@code >>}）の後ろで {@code 123} を
 * 入力し, {@code IntValue} の実行で「結果」が 123 になり, 続いて次のプロンプトの後ろで
 *  {@code + 45} を入力し, {@code IntAdd} の実行で「結果」が 123 + 45 の結果である
 * 168 になったものである. 
 */
class IntCalc {
    /**
     * 電卓を作って実行する. 
     */
    public static void main(String [] args) throws Exception {
        // 定義した様々な「コマンド」からなるリストを作る
        ArrayList<Command<BigInteger>> comms = new ArrayList<Command<BigInteger>>();
        comms.add(new EmptyCommand<BigInteger>());
        comms.add(new IntValue());
        comms.add(new IntAdd());
        comms.add(new IntSub());
        comms.add(new IntMulDiv());
        comms.add(new IntNeg());
        // 標準入力から読み込む BufferedReader を作る
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        // 以上の二つを与えて, 新たな電卓のインスタンスを生成
        Calculator<BigInteger> c = new Calculator<BigInteger>(br, comms);
        // 電卓の実行
        c.run(BigInteger.ZERO);
    }
}
