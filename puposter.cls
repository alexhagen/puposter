\NeedsTeXFormat{LaTeX2e}[1998/06/01]
\ProvidesClass{puposter}[2016/06/22 class for purdue poster design]

%% Load base
\LoadClass{article}

\RequirePackage{fontspec}
\RequirePackage{tikz} % to draw
\RequirePackage{anyfontsize} % to use the fonts
\RequirePackage{lipsum}
\RequirePackage{pst-all}
\RequirePackage{pst-light3d}
\usetikzlibrary{backgrounds}

\pgfdeclarelayer{background}
\pgfdeclarelayer{lines}
\pgfdeclarelayer{foreground}
\pgfsetlayers{background,lines,main,foreground}

% Change the paper to 3' x 4'
\RequirePackage[paperheight=36in,
                paperwidth=48in,
                margin=1in]{geometry}

\newfontfamily\champion[Path = fonts/]{champion-htf-featherweight.otf}
\newfontfamily\myriad[Path = fonts/,%
                      ItalicFont=MyriadPro-CondIt.otf,%
                      BoldFont=MyriadPro-BoldCond.otf,%
                      BoldItalicFont=MyriadPro-BoldCondIt.otf%
                     ]{MyriadPro-Cond.otf}
\newfontfamily\chaparral[Path=fonts/,%
                         ItalicFont=ChaparralPro-Italic.otf,%
                         BoldFont=ChaparralPro-Bold.otf,%
                         BoldItalicFont=ChaparralPro-BoldIt.otf%
                        ]{ChaparralPro-Regular.otf}
\setmainfont[Path=fonts/,%
             ItalicFont=MyriadPro-CondIt.otf,%
             BoldFont=MyriadPro-BoldCond.otf,%
             BoldItalicFont=MyriadPro-BoldCondIt.otf%
            ]{MyriadPro-Cond.otf}
\setsansfont[Path=fonts/,%
             ItalicFont=MyriadPro-CondIt.otf,%
             BoldFont=MyriadPro-BoldCond.otf,%
             BoldItalicFont=MyriadPro-BoldCondIt.otf%
            ]{MyriadPro-Cond.otf}

%\DeclareMathSizes{10}{14}{10}{6}
%\DeclareMathSizes{12}{14}{10}{6}
%\DeclareMathSizes{14.4}{20}{14}{10}

\definecolor{grey20}{RGB}{209,211,212}
\definecolor{grey40}{RGB}{167,169,172}
\definecolor{grey60}{RGB}{116,108,102}
\definecolor{oldgold}{RGB}{163,121,44}
\definecolor{newgold}{RGB}{227,174,36}

%\usetikzlibrary{positioning}
\usetikzlibrary{calc}
\usetikzlibrary{fadings}
\usetikzlibrary{patterns}

\newcommand{\conference}[1]{\def\@conference{#1}}
\newcommand{\affiliation}[1]{\def\@affiliation{#1}}
\newcommand{\tss}[1]{\textsuperscript{\myriad\mdseries{#1}}}

\input{pupatterns.tex}

\pagestyle{empty}

\renewcommand\maketitle{
    \begin{tikzpicture}[remember picture, overlay]
        \begin{pgfonlayer}{lines}
            \fill[newgold]
                  ($(current page.north west) + (1.5in, -1.5in)$)
                  rectangle
                  ($(current page.north east) + (-0.5in, -7.0in)$);
            \fill[hatchcolor=white, pattern=custom north west lines,
                  hatchspread=50pt, hatchthickness=20pt]
                  ($(current page.north west) + (1.0in, -1.0in)$)
                  rectangle
                  ($(current page.north east) + (-0.5in, -7.5in)$);

            \fill[grey40]
                ($(current page.north west) + (1in, -1in)$) --
                ($(current page.north east) + (-1in, -1in)$) --
                ($(current page.north east) + (-1in, -6.5in)$) --
                ($(current page.north west) + (1in, -6.5in)$) --
                cycle;
        \end{pgfonlayer}
        \begin{pgfonlayer}{main}
            \node[
              align=left,
              font=\linespread{0.7}\champion\color{newgold}\fontsize{140}{156}\mdseries,
              anchor=north west,
              minimum width=\paperwidth - 2.75in,
              text width=\paperwidth - 2.75in, inner sep=0, outer sep=0
              ]
              (title) at ($(current page.north west) + (1.5in, -1.75in)$)
              {\@title};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{56}{64}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (authors) at ($(title.east) + (-0.5in, -.75in)$)
            {\@author};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{56}{64}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (conference) at ($(authors.south east) + (0in, -0.25in)$)
            {\@conference};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{36}{42}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (affiliation) at ($(current page.north east) + (-0.5in, -7.125in)$)
            {\@affiliation};
        \end{pgfonlayer}
    \end{tikzpicture}
}

\AtBeginDocument{
\begin{tikzpicture}[remember picture,overlay]
    \begin{pgfonlayer}{background}
        \node[scale=1.5,
              inner sep=0,
              outer sep=0,
              anchor=south east]
             at (current page.south east)
             {\includegraphics{css/background_bottom_right.pdf}};
        \node[scale=1.5,
              inner sep=0,
              outer sep=0,
              anchor=north west]
             at (current page.north west)
             {\includegraphics{css/background_top_left.pdf}};
    \end{pgfonlayer}
\end{tikzpicture}
}
