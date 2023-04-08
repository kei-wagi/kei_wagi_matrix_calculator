/*
 * 電卓のコア部分の実装.  
 * ドキュメントの生成はこのコマンドで：
 *
 *   javadoc -private -d doc Calculator.java *Calc.java
 * 
 */

import java.util.*;
import java.io.*;
import java.math.*;

/**
 * テキストベースの電卓クラス. 
 * 電卓は, 現在の「結果」をひとつ保持しており, ユーザがターミナルに入力した「ブロック」に応じて演算を実行し, 
 * 「結果」を更新する. <br />
 * 電卓で実行する演算（「コマンド」）は, 
 * {@code Command} インターフェースを実装したオブジェクトが提供する. <br />
 * 電卓は, 「コマンド」のリスト（コマンドリスト）を保持しており, 
 * ユーザの入力した「ブロック」に対して, 
 * 実行可能な「コマンド」を実行して現在の「結果」を更新する. <br />
 * 「ブロック」とは, トークンの並んだ一行
 * <p><blockquote><pre>
 * token_1 token_2 ... token_k
 * </pre></blockquote><p>
 * か, 行末にコロン ':' が置かれた一行からタブ（下では TAB で表現）でインデントされた行が続く限り
 * <p><blockquote><pre>
 * token_1 token_2 ... token_k :
 *  TAB line_1 
 *  TAB line_2
 *     ...
 *  TAB line_n 
 * </pre></blockquote><p>
 * のいずれかのことである. なお, 後者の場合には, 「ブロック」の終端として空行が入力されるとする. <br />
 * 例えば, 通常の電卓のように「現在の値に 10 を足す」を行うなら
 * <p><blockquote><pre>
 * + 10
 * </pre></blockquote><p>
 * のような「ブロック」をユーザは入力することになるだろう. 
 * また, 例えば行列電卓で 2x2 行列を足すのであれば, 
 * <p><blockquote><pre>
 * add :
 *  TAB  1 2
 *  TAB  3 4
 * </pre></blockquote><p>
 * のように入力させることになるだろう. 
 * @param Result 電卓の「結果」の型. 
 */
class Calculator<Result> {
    /**
     * ユーザ入力を読み込むための {@code BufferedReader}. 
     * 標準入力なのか実際のファイルなのかは, この電卓としては気にしない. 
     */
    BufferedReader br;
    /**
     * 「コマンド」のリスト. 
     * この電卓自身は演算を持っておらず, 外部から「コマンド」のリストを受け取ることで実行可能な演算が決まる. 
     */
    List<Command<Result>> comms;
    /**
     * 与えられた {@code BufferedReader} から入力を読み込み, 
     * 与えられた「コマンド」のリストにある演算を実行する電卓を作るコンストラクタ. 
     * @param br ここから入力を行単位で読み込む. 
     * @param comms このリストにある「コマンド」を電卓の演算とする. 
     */
    Calculator(BufferedReader br, List<Command<Result>> comms) {
        this.br = br;
        this.comms = comms;
    }
    /**
     * 電卓のメインループ.
     * ユーザの入力する「ブロック」をひとつずつ読み込み, 
     * それをコマンドリストの各「コマンド」に順に実行を問い合わせ, 
     * 「コマンド」が「結果」を返したら次の「ブロック」の処理に向かう. <br />
     * 入力が尽きたらループを終了してその時点の「結果」を返す. 
     * @param res 電卓の初期値とする「結果」. 
     * @return 電卓の最終的な「結果」. 
     */
    Result run(Result res) {
        showCurrentResult(res);   // とりあえず最初に現在の「結果」を表示
        for(;;) {
            List<String> block = getNextBlock();   // ユーザの入力した「ブロック」を取得
            if(block == null) break;               // 入力が尽きたらループ終了
            boolean run = false;                   // 「コマンド」が実行できたか？
            String [] ts = tokenize(block.get(0)); // 1行目をトークンに分解
            for(Command<Result> c : comms) {       // 各「コマンド」について
                Result r = c.tryExec(ts, block, res);  // 実行を問い合わせる
                if(r != null) {                        // null 以外の値がきた → 実行できた
                    res = r;                           // その値を「結果」に保存し、表示
                    showCurrentResult(res);
                    run = true;
                    break;
                }
            }
            // ひとつの「コマンド」も実行できなかった → エラー表示して続行
            if(!run) {
                System.err.println("Unknown command: \"" + block.get(0) + "\"");
            }
        }
        return res;
    }
    /**
     * 文字列をトークンに分割する.
     * トークンの切れ目になるのは, 1文字以上の空白文字か, 記号の前後かである.
     * よって, {@code 10 + x*y} という文字列は, 
     * {@code 10} と {@code +} と {@code x} と {@code *} と {@code y} という5個のトークンに切られる. 
     * @param line 分割対象の文字列.
     * @return 分割されたトークンが順に並んだ配列.
     */
    String [] tokenize(String line) {
        // 最初に記号の前後に空白を入れてから, 空白文字でぶった切る
        return line.replaceAll("(\\W)"," $1 ").replaceAll("^\\s+","").split("\\s+");
    }
    /**
     * 与えられた「結果」を標準出力へ表示する. 
     * @param res 表示したい「結果」.
     */
    void showCurrentResult(Result res) {
        System.out.println(res);
    }
    /**
     * 入力から「ブロック」をひとつ切り出す.
     * 「ブロック」は, 行（文字列）のリストとして返される.
     * 「ブロック」の終端を表すための空行の入力はリストに含まれない. 
     * また，2行目以降の先頭の TAB も削除されている（先頭に複数の TAB が連続していても，最初の TAB ひとつのみを削除）．
     * @return 「ブロック」を構成する行（文字列）のリスト. 入力が尽きている場合には {@code null}. 
     */
    List<String> getNextBlock() {
        try {
            List<String> block = new ArrayList<String>();
            String line = readNextLine();
            if(line == null) {   // 入力が尽きた場合は null
                return null;
            }
            // 末尾が : なら複数行バージョンの「ブロック」
            if(line.length() > 0 && line.charAt(line.length()-1) == ':') { 
                block.add(line.substring(0, line.length() - 1)); // 末尾の : は削っとく
                for(;;) {
                    line = readNextBlockLine();
                    // 次の行の先頭がタブなら, リストに追加. そうでなければ読み込み終了. 
                    if(line.length() == 0 || line.charAt(0) != '\t') {
                        if(line.length() > 0) {
                            System.err.println("Warn: ignoring extra line: " + line);
                        }
                        break;
                    }
                    block.add(line.substring(1)); // タブ文字は削っておく
                }
            } else { // 1行のみの「ブロック」の場合
                block.add(line);
            }
            return block;
        } catch(IOException e) { // なにか例外が起きたら入力が終わったとする
            return null;
        }
    }
    /**
     * ユーザに入力を求めるプロンプトの文字列表示用. 
     * @param str プロンプトとして表示する文字列.
     */
    void showPrompt(String str) {
        System.out.print(str);
    }
    /**
     * 入力から「ブロック」の先頭1行を読み込む. 
     * プロンプトを表示して {@code br} から1行読み込む.
     * 読み込んだ行の先頭に TAB があったら警告して空行を返す. 
     * @return 読み込まれた一行. 
     */
    String readNextLine() throws IOException {
        showPrompt(">> ");
        String line = br.readLine();
        if(line != null && line.length() > 0 && line.charAt(0) == '\t') {
            System.err.println("Unexpected TAB is found: " + line);
            return "";
        }
        return line;
    }
    /**
     * 入力から「ブロック」の2行目以降を読み込む. 
     * @return 読み込まれた一行. 
     */
    String readNextBlockLine() throws IOException {
        showPrompt(".. ");
        return br.readLine();
    }
}

/**
 * 電卓の「コマンド」のインターフェース. 
 * 「コマンド」は, 電卓の現在の「結果」を受け取り新たな「結果」を作るもの. 
 * @param Result 電卓の「結果」の型. 
 */
interface Command<Result> {
    /**
     * 「コマンド」の実行を行うメソッド. 
     * 与えられた「ブロック」がこの「コマンド」の実行に対応しているなら, 
     * 受け取った現在の「結果」と「ブロック」から新たな「結果」を計算して返す. 
     * 「ブロック」がこの「コマンド」の実行でないならば, {@code null} を返す. 
     * @param ts 「ブロック」の最初の行をトークンに分解して並べた配列. 
     * @param block 「ブロック」全体
     * @param res 電卓の現在の「結果」
     * @return 「ブロック」がこの「コマンド」のものであれば, 新たな「結果」を返す. 
               そうでない場合, {@code null}. 
     */
    Result tryExec(final String [] ts, final List<String> block, final Result res);
}


/**
 * 空行が入力されたときに何もしないための「コマンド」. 
 * @param Result 電卓の「結果」の型. 
 */
class EmptyCommand<Result> implements Command<Result> {
    /**
     * 「ブロック」が1行のみで, かつ, その行が空（空白文字のみ）のときに実行され, 現在の「結果」をそのまま返す. 
     * @return 「ブロック」が条件に合っていれば, 与えられた「結果」そのもの. そうでない場合, {@code null}.
     */
    public Result tryExec(final String [] ts, final List<String> block, final Result res) {
        if(block.size() == 1 && ts.length == 1 && ts[0].equals("")) return res;
        return null;
    }
}



