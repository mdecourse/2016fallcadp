<?xml version="1.0" encoding="utf-8"?>
<!-- Created by Leo: http://leoeditor.com/leo_toc.html -->
<?xml-stylesheet ekr_test ?>
<leo_file xmlns:leo="http://leoeditor.com/namespaces/leo-python-editor/1.1" >
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5" body_secondary_ratio="0.5">
	<global_window_position top="50" left="50" height="500" width="700"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="ibm.20161021125730.1"><vh>@settings</vh>
<v t="ibm.20161021125730.2"><vh>@data qt-gui-plugin-style-sheet</vh></v>
<v t="ibm.20161021125730.3"><vh>@string initial_split_orientation = horizontal</vh></v>
</v>
<v t="ibm.20161021125710.2"><vh>有關 2016fallcadp</vh></v>
<v t="ibm.20161021125827.1"><vh>刪除 gh-page 分支</vh></v>
<v t="ibm.20161021130015.1" a="E"><vh>課程議題</vh>
<v t="ibm.20161021130112.1" a="E"><vh>期中自評連結</vh>
<v t="ibm.20161021130311.1"><vh>2a分組名單</vh></v>
<v t="ibm.20161021130236.1"><vh>2b分組名單</vh></v>
<v t="ibm.20161021130401.1"><vh>組員如何互評</vh></v>
</v>
</v>
</vnodes>
<tnodes>
<t tx="ibm.20161021125710.2">2016 Fall CADP 電腦輔助設計實習課程

倉儲位置: https://github.com/mdecourse/2016fallcadp

每一協同成員都應該要在 users 目錄下有一個 Leo 檔案

</t>
<t tx="ibm.20161021125730.1"></t>
<t tx="ibm.20161021125730.2">QTreeWidget {
    /* These apply to the selected item, but not to editing items.*/
    background-color: #ffffec; /* Leo's traditional tree color */
    selection-color: black; /* was white */
    selection-background-color: lightgrey;
    /* font-family: SansSerif; */
    /*font-family: DejaVu Sans Mono;*/
    font-family:YaHei Mono;
    /* 標題字型大小設定 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
 }

/* Headline edit widgets */
QTreeWidget QLineEdit {
    background-color: cornsilk;
    selection-color: white;
    selection-background-color: blue;
    /*font-family: DejaVu Sans Mono;*/    
    font-family:YaHei Mono;
    /* 沒有特別對應字型大小 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
}

/* The log panes */
QTextEdit {
    background-color: #f2fdff;
    selection-color: red;
    selection-background-color: blue;
    /* font-family: Courier New; */
    font-family:YaHei Mono;
    /* log font 大小 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
}

/* The body pane */
QTextEdit#richTextEdit {
    background-color: #fdf5f5; /* A kind of pink. */
    selection-color: white;
    selection-background-color: red;
    /*font-family: DejaVu Sans Mono;*/
    /* font-family: Courier New; */
    font-family:YaHei Mono;
    /* 內文字型大小 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
}

QLabel {
    font-family:YaHei Mono;'CherryPy', 'pytz', 'mako', 'beautifulsoup4', 'pymysql', 'peewee'
    /* 下方的 Minibuffer 標題字型大小 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
}

/* Editor labels */
QLineEdit#editorLabel {
    background-color: #ffffec;
    font-family:YaHei Mono;
    /* 沒有直接對應字型大小 */
    font-size: 20px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
    border: 2px;
    margin: 2px;
}</t>
<t tx="ibm.20161021125730.3">horizontal: body pane to the left
vertical: body pane on the botton</t>
<t tx="ibm.20161021125827.1">git branch -d local_branch</t>
<t tx="ibm.20161021130015.1">電腦使用
網路使用
github 使用
Vimeo 使用
Solvespace 使用
OnShape 使用
</t>
<t tx="ibm.20161021130112.1">已知各班分組名單, 必須設法建立班上網頁檢視與各組網頁檢視與各學員網頁檢視連結</t>
<t tx="ibm.20161021130236.1">第1組:40123252,40223128,40223130,40223215,40223233,
第2組:40223205,40223206,40423225,40423231,40423240,40423250,
第3組:40423201,40423202,40423203,40423211,40439221,42423223,
第4組:40423204,40423221,40423228,40423243,40423253,40423254,
第5組:40423205,40423227,40423235,40423236,40423237,40423238,
第6組:40423206,40423214,40423216,40423224,40423226,40423229,
第7組:40423207,40423210,40423219,40423220,40423241,40423246,
第8組:40423208,40423213,40423232,40423233,
第9組:40423209,40423212,40423217,40423218,40423239,40423255,
第10組:40423222,40423245,40423248,40423251,
第11組:40423234,40423242,40423244,40423247,40423249,40423252,</t>
<t tx="ibm.20161021130311.1">第1組:40223103,40223109,40423111,40443111,
第2組:40423101,40423103,40423104,40423117,40423136,40423144,
第3組:40423105,40423115,40423119,40423128,40423139,40423149,
第4組:40423106,40423107,40423108,40423121,40423125,40423140,
第5組:40423109,40423110,40423113,40423130,40423131,40423138,
第6組:40423112,40423116,40423118,40423127,40423137,40423152,
第7組:40423114, 40423120,40423122,40423124,
第8組:40423123,40423129,40423132,40423133,40423155,40423157,
第9組:40423126,40423135,40423150,40423154,
第10組:40423141,40423143,40423145,40423146,40423147,40423148,</t>
<t tx="ibm.20161021130401.1">各組員除自評外, 可以進行互評, 並各自簡要說明評分理由</t>
</tnodes>
</leo_file>
