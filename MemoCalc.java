/*
 * メモ機能付きの整数電卓. 
 * コンパイル & 実行：
 * javac Calculator.java IntCalc.java MemoCalc.java
 * java MemoCalc
 */

import java.util.*;
import java.io.*;
import java.math.*;

/**
 * 変数に「結果」を記憶しておくメモリを表すクラス. 
 * なお, これ自身が「保存している変数の一覧を表示する」という「コマンド」になっている. 
 * <p><blockquote><pre>{@code
 * show
 * }</pre></blockquote><p>
 * という1トークンからなる1行「ブロック」を受け付ける. 
 * @param Result 変数に入れる「結果」の型
 */
class Memory<Result> implements Command<Result> {
    /**
     * 変数名（String）から「結果」への連想配列. これに変数の情報を保存する. 
     */
    HashMap<String, Result> mem;
    /**
     * 何の変数も存在していないメモリを作るコンストラクタ. 
     */
    Memory() {
        mem = new HashMap<String, Result>();
    }
    /**
     * 変数に保存されている「結果」を返す. 
     * 実際のところは連想配列に問い合わせるだけ. 
     * @param var 変数名
     * @return その変数に保存された「結果」. なければ {@code null}
     */
    public Result get(String var) {
        return mem.get(var);
    }
    /**
     * 変数に与えられた「結果」を保存する. 
     * 実際のところは連想配列にエントリを追加するだけ. 
     * @param var 変数名
     * @param val その変数に保存する「結果」
     */
    public void put(String var, Result val) {
        mem.put(var, val);
    }
    /**
     * 保存されている変数の一覧を標準出力に表示する. 
     */
    public Result tryExec(final String [] ts, final List<String> block, final Result res) {
        if(block.size() != 1) return null;
        if(ts.length == 1 && "show".equals(ts[0])) {
            for(String var : mem.keySet()) {
                Result v = mem.get(var);
                String sv = v.toString();
                if(sv.indexOf('\n') >= 0) {
                    sv = sv.replaceAll("^|\\n", "\n ");
                }
                System.out.println(var + " = " + sv);
            }
            return res;
        }
        return null;
    }
}

/**
 * 変数の値を必要とする「コマンド」のためのベースクラス. 
 * いまのところほとんど意味はない. 
 * 変数の値を必要とする「コマンド」だけに対してなにか統一的な処理をしたいときに役立つだろう. 
 * @param Result 「結果」の型（＝変数の値の型）
 */
abstract class CommandWithMemory<Result> implements Command<Result> {
    /**
     * 変数の値を保持しておくオブジェクト. 
     */
    Memory<Result> mem;
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    CommandWithMemory(Memory<Result> mem) {
        this.mem = mem;
    }
}

/**
 * 変数に現在の「結果」を保存し, また, 変数に保存した値を現在の「結果」に読み出す「コマンド」. 
 * <p><blockquote><pre>{@code
 * store var
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, 現在の「結果」をその変数 {@code var} に保存する. 
 * また, 
 * <p><blockquote><pre>{@code
 * load var
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, その変数 {@code var} に保存された値を現在の「結果」にする. 
 */
class LoadStore<Result> extends CommandWithMemory<Result> {
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    LoadStore(Memory<Result> mem) {
        super(mem); // 親コンストラクタをそのまま呼ぶだけ
    }
    /**
     * 変数への保存と読み出しを実行する. 
     * 読み出し時に変数が見つからない場合には, {@code UnknownVariableException} 例外を投げる. 
     */
    public Result tryExec(final String [] ts, final List<String> block, final Result res) {
        if(block.size() != 1) return null;
        if(ts.length != 2) return null;
        if("load".equals(ts[0])) {
            Result m = mem.get(ts[1]);
            if(m == null) throw new UnknownVariableException(ts[1]);
            return m;
        }
        if("store".equals(ts[0])) {
            mem.put(ts[1], res);
            return res;
        }
        return null;
    }
}

/**
 * 変数が見つからなかったときに使うための例外. 
 * 面倒なので {@code RuntimeException} としておく. 
 */
class UnknownVariableException extends RuntimeException {
    /**
     * 見つからなかった変数名を受け取るコンストラクタ. 
     * @param var 見つからなかった変数. 
     */
    UnknownVariableException(String var) {
        super("Unknown variable: " + var);
    }
}

/**
 * 加減乗除のための「コマンド」. 
 * <p><blockquote><pre>{@code
 * op n
 * }</pre></blockquote><p>
 * ないし
 * <p><blockquote><pre>{@code
 * op var
 * }</pre></blockquote><p>
 * という 1行の「ブロック」を受け付けて, その値{@code n}ないし変数{@code var}の値を
 * 現在の「結果」に加減乗除した「結果」を返す. 
 * {@code op} は {@code +} か {@code -} か {@code *} か {@code /} であり, 
 * 対応する演算が行われる. 
 */
class IntArithWithMemory extends CommandWithMemory<BigInteger> {
    /**
     * 変数の情報を保持する {@code Memory} オブジェクトを受け取るコンストラクタ. 
     * @param mem 変数の情報を保持するオブジェクト. 
     */
    IntArithWithMemory(Memory<BigInteger> mem) {
        super(mem); // 親クラスのコンストラクタそのままよぶ
    }
    /**
     * トークンをひとつうけとり, 変数の値の読み出しか 10進数としての解釈をこの順で試して値を返す. 
     * なお, 手抜き実装にしているので 10進数として解釈できないと例外が投げられる. 
     * @param token トークン. 変数名であるか, 10進数であることが期待される. 
     * @return トークンが保存された変数名であった場合にはその変数の値を, 
               そうでなく, 10進数であった場合にはその数値を返す. それ以外の場合には {@code null}.
     */
    BigInteger eval(String token) {
        BigInteger m = mem.get(token); // 変数の値を Memory に問い合わせる
        if(m == null) m = new BigInteger(token); // 変数の値がなければ 10進数として解釈
        return m;
    }
    /**
     * 1行「ブロック」の最初のトークンに指定された演算子での演算を, 2つ目のトークン（変数か値）で実行する. 
     * 
     */
    public BigInteger tryExec(final String [] ts, final List<String> block, final BigInteger res) {
        if(block.size() != 1) return null;
        if(ts.length == 2 && "+".equals(ts[0])) {
            return res.add(eval(ts[1]));
        } else if(ts.length == 2 && "-".equals(ts[0])) {
            return res.subtract(eval(ts[1]));
        } else if(ts.length == 2 && "*".equals(ts[0])) {
            return res.multiply(eval(ts[1]));
        } else if(ts.length == 2 && "/".equals(ts[0])) {
            return res.divide(eval(ts[1]));
        } else {
            return null;
        }
    }
}

/**
 * 行列電卓を作成して動作させるクラス. 
 * 例えば, ターミナルで次のような実行ができる. 
 * <p><blockquote><pre>{@code
 * $ javac Calculator.java IntCalc.java MemoCalc.java
 * $ java MemoCalc
 * 0
 * >> + 10
 * 10
 * >> store x
 * 10
 * >> + x
 * 20
 * >> store y
 * 20
 * >> show
 * x = 10
 * y = 20
 * 20
 * >> load x
 * 10
 * >>
 * }</pre></blockquote><p>
 * これは, 最初に「結果」が 0 で電卓が動き始め, 
 * まずは最初のプロンプトの後ろで {@code + 10} と（1行の「ブロック」を）入力し, 
 * その「ブロック」に対して {@code IntArithWithMemory} が実行されて「結果」が 10 になり, 
 * 続いて {@code store x} と入力し, 
 * {@code LoadStore} が動作してその値（10）を変数 x に保存し, 
 * 続いて {@code + x} と入力し, で現在の「結果」に変数 x の値（10）が加算されて「結果」が 20 になり, 
 * さらに {@code store x} と入力し, {@code LoadStore} が動作してその値（20）を変数 y に保存し, 
 * そこで {@code show} と入力して {@code Memory} が動作して x と y の値が表示され, 
 * 最後に {@code load x} と入力し, 現在の「結果」が変数 x の値である 10 になった. 
 */
class MemoCalc {
    /**
     * メモリ機能付き電卓の生成と実行. 
     */
    public static void main(String [] args) throws Exception {
        // 変数の値を覚えておくメモリオブジェクトの生成
        Memory<BigInteger> mem = new Memory<BigInteger>();
        // コマンドリストを生成. 一部の「コマンド」は, メモリオブジェクトを共有する
        ArrayList<Command<BigInteger>> comms = new ArrayList<Command<BigInteger>>();
        comms.add(new EmptyCommand<BigInteger>());
        comms.add(new IntValue());
        comms.add(new IntNeg());
        comms.add(new IntArithWithMemory(mem));
        comms.add(new LoadStore<BigInteger>(mem));
        comms.add(mem);
        // 入力は標準入力から
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        // 電卓オブジェクトの生成と実行
        Calculator<BigInteger> c = new Calculator<BigInteger>(br, comms);
        c.run(BigInteger.ZERO);
    }
}

