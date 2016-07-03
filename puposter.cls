\NeedsTeXFormat{LaTeX2e}[1998/06/01]
\ProvidesClass{puposter}[2016/06/22 class for purdue poster design]

%% Load base
\LoadClass[12pt]{article}

\RequirePackage{fontspec}
\RequirePackage{calc}
\RequirePackage{tikz}
\RequirePackage{tikzpagenodes}
\RequirePackage{anyfontsize}
\RequirePackage{esint}
\RequirePackage{lipsum}
\RequirePackage{pst-all}
\RequirePackage{pst-light3d}
\RequirePackage{totcount}
\RequirePackage{environ}
\usetikzlibrary{shadows}
\tikzset{
    text shadow/.code args={[#1]#2at#3(#4)#5}{
        \pgfkeysalso{/tikz/.cd,#1}%
        \foreach \angle in {0,5,...,359}{
                \node[#1,text=white] at ([shift={(\angle:5pt)}] #4){#5};
        }
    }
}

\usetikzlibrary{backgrounds}
\pgfdeclarelayer{background}
\pgfdeclarelayer{lines}
\pgfdeclarelayer{foreground}
\pgfsetlayers{background,lines,main,foreground}

% Change the paper to 3' x 4'
\RequirePackage[paperheight=36in,
                paperwidth=48in,
                margin=1in,
                top=8in,
                bottom=3.5in,
                text={46in, 25in}]{geometry}

\newfontfamily\champion[Path = fonts/,%
                        Scale = 2.0%
                       ]{champion-htf-featherweight.otf}
\newfontfamily\myriad[Path = fonts/,%
                      ItalicFont=MyriadPro-CondIt.otf,%
                      BoldFont=MyriadPro-BoldCond.otf,%
                      BoldItalicFont=MyriadPro-BoldCondIt.otf,%
                      Scale = 2.0%
                     ]{MyriadPro-Cond.otf}
\newfontfamily\chaparral[Path=fonts/,%
                         ItalicFont=ChaparralPro-Italic.otf,%
                         BoldFont=ChaparralPro-Bold.otf,%
                         BoldItalicFont=ChaparralPro-BoldIt.otf,%
                         Scale = 2.0%
                        ]{ChaparralPro-Regular.otf}
\setmainfont[Path=fonts/,%
             ItalicFont=MyriadPro-CondIt.otf,%
             BoldFont=MyriadPro-BoldCond.otf,%
             BoldItalicFont=MyriadPro-BoldCondIt.otf,%
             Scale = 2.0%
            ]{MyriadPro-Cond.otf}
\setsansfont[Path=fonts/,%
             ItalicFont=MyriadPro-CondIt.otf,%
             BoldFont=MyriadPro-BoldCond.otf,%
             BoldItalicFont=MyriadPro-BoldCondIt.otf,%
             Scale = 2.0%
            ]{MyriadPro-Cond.otf}

\DeclareMathSizes{9.8}{24}{12}{12}
\DeclareMathSizes{10.0}{24}{12}{12}
\DeclareMathSizes{10.95}{24}{12}{12}   % For size 10 text
\DeclareMathSizes{11}{24}{18}{14}      % For size 11 text
\DeclareMathSizes{12}{30}{24}{16}     % For size 12 text

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
                  ($(current page.north east) + (-0.5in, -7.75in)$);

            \fill[grey40]
                ($(current page.north west) + (1in, -1in)$) --
                ($(current page.north east) + (-1in, -1in)$) --
                ($(current page.north east) + (-1in, -6.5in)$) --
                ($(current page.north west) + (1in, -6.5in)$) --
                cycle;
        \end{pgfonlayer}%
        \begin{pgfonlayer}{main}
            \node[
              align=left,
              font=\linespread{0.7}\champion\color{newgold}\fontsize{70}{78}\mdseries,
              anchor=north west,
              minimum width=\paperwidth - 2.75in,
              text width=\paperwidth - 2.75in, inner sep=0, outer sep=0
              ]
              (title) at ($(current page.north west) + (1.5in, -1.75in)$)
              {\@title};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{28}{32}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (authors) at ($(title.east) + (-0.5in, -.75in)$)
            {\@author};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{28}{32}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (conference) at ($(authors.south east) + (0in, -0.25in)$)
            {\@conference};
            \node[
              align=right,
              font=\myriad\color{black}\fontsize{18}{22}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (affiliation) at ($(current page.north east) + (-0.5in, -7.125in)$)
            {\@affiliation};
        \end{pgfonlayer}
    \end{tikzpicture}%
    \noindent%
}

\AtBeginDocument{
    \setlength\parindent{0pt}%
    \SetMathAlphabet\mathrm{normal}{OT1}{myriadmath}{m}{n}%
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
        \end{pgfonlayer}%
        \begin{pgfonlayer}{lines}
            \draw[draw=grey60] ($(current page.south west) + (1.5in, 3.0in)$) --
                               ($(current page.south east) + (-1.5in, 3.0in)$);
        \end{pgfonlayer}%
    \end{tikzpicture}%
    \def\incell{0}%
}

\newcommand{\outertitle}[1]{
  \begin{tikzpicture}
      \node[inner sep=0.0in, outer sep=0.0in, text width=\textwidth,
            fill=grey40, align=center, anchor=west, minimum height=1.25in,
            text depth=0.2in]
      {
          \color{newgold}\fontsize{36}{42}\champion{#1}
      };
  \end{tikzpicture}%
}

\newcommand{\innertitle}[2]{
  \begin{tikzpicture}
      \node[fill=newgold, minimum width=\textwidth,
            minimum height=1.5in] (goldrect) at (0, 0) {};
      \node[hatchcolor=white, pattern=custom north west lines,
            hatchspread=20pt, hatchthickness=7pt, hatchshift=5pt,
            minimum width=\textwidth + 1in, minimum height=2in,
            anchor=center]
            (linesrect) at ($(goldrect.center) + (0in,0in)$) {};
      \node[anchor=center, text width=\textwidth, outer sep=0in, inner sep=0.125in,
            text shadow={[align=center] at (goldrect.center) {\Huge \champion #1\\ \LARGE \champion #2}}]
      at (goldrect.center)
            {\Huge \champion \color{newgold} #1\\ \LARGE \champion #2};
  \end{tikzpicture}%
}

\newlength{\cellheight}
\newlength{\cellwidth}
\NewEnviron{cell}[4]
{
  \noindent%
  \setlength{\cellheight}{2.0in * \real{#2}}%
  \setlength{\cellwidth}{3.5in * \real{#1}}%
  \begin{minipage}[t][\cellheight][t]{\cellwidth}
    \noindent%
    % Check if we're already inside a cell, if not, we make the outer title
    \ifnum\incell=1
      \ifx&#3&%
        % #1 is empty
      \else
        \innertitle{#3}{#4}%
        \vspace{0.06125in}
      \fi
    \else
      \outertitle{#3}%
      \vspace{0.125in}
    \fi
    \def\incell{1}%
    \BODY%
  \end{minipage}%
}
