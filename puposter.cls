\NeedsTeXFormat{LaTeX2e}[1998/06/01]
\ProvidesClass{puposter}[2016/06/22 class for purdue poster design]

%% Load base
\LoadClass[12pt]{article}

\RequirePackage[no-math]{fontspec}
\RequirePackage{calc}
\RequirePackage{array}
\RequirePackage[table]{xcolor}
\RequirePackage{tikz}
\RequirePackage{fp}
\RequirePackage{tikzpagenodes}
\RequirePackage{svg}
\RequirePackage{anyfontsize}
\RequirePackage{adjustbox}
\RequirePackage{ifthen}
\RequirePackage{lipsum}
\RequirePackage{float}
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

\RequirePackage{mathastext}

%\DeclareMathSizes{9.8}{24}{12}{12}
%\DeclareMathSizes{10.0}{24}{12}{12}
%\DeclareMathSizes{10.95}{24}{12}{12}   % For size 10 text
%\DeclareMathSizes{11}{24}{14}{10}      % For size 11 text
%\DeclareMathSizes{12}{26}{18}{12}     % For size 12 text

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
            \fill[grey40]
                  ($(current page.north west) + (1.5in, -1.5in)$)
                  rectangle
                  ($(current page.north east) + (-0.5in, -7.0in)$);
            \fill[hatchcolor=white, pattern=custom north west lines,
                  hatchspread=50pt, hatchthickness=20pt]
                  ($(current page.north west) + (1.0in, -1.0in)$)
                  rectangle
                  ($(current page.north east) + (-0.5in, -7.75in)$);

            \fill[newgold]
                ($(current page.north west) + (1in, -1in)$) --
                ($(current page.north east) + (-1in, -1in)$) --
                ($(current page.north east) + (-1in, -6.5in)$) --
                ($(current page.north west) + (1in, -6.5in)$) --
                cycle;
        \end{pgfonlayer}%
        \begin{pgfonlayer}{main}
            \node[
              align=left,
              font=\linespread{0.7}\champion\color{black}\fontsize{70}{78}\mdseries,
              anchor=north west,
              minimum width=\paperwidth - 2.75in,
              text width=\paperwidth - 2.75in, inner sep=0, outer sep=0
              ]
              (title) at ($(current page.north west) + (1.5in, -1.75in)$)
              {\@title};
            \node[
              align=right,
              font=\myriad\color{grey60}\fontsize{28}{32}\mdseries,
              anchor=north east,
              inner sep=0,
              outer sep=0
            ]
            (authors) at ($(title.east) + (-0.5in, -.75in)$)
            {\@author};
            \node[
              align=right,
              font=\myriad\color{grey60}\fontsize{28}{32}\mdseries,
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
    \def\numheadings{0}
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
  \begin{tikzpicture}[remember picture, overlay]
  \begin{pgfonlayer}{lines}
      \node[inner sep=0.0in, outer sep=0.0in, text width=\textwidth,
            fill=newgold, align=center, anchor=north west, minimum height=1.25in,
            text depth=0.2in]
      {
          \color{black}\fontsize{36}{42}\champion{#1}
      };
  \end{pgfonlayer}
  \end{tikzpicture}%
}

\newcommand{\innertitle}[2]{
  \begin{tikzpicture}[remember picture, overlay]
    \begin{pgfonlayer}{background}
      \node[fill=grey40, minimum width=\textwidth,
            minimum height=1.5in, anchor=north west] (goldrect) at (0.125in, 0) {};
      \node[hatchcolor=white, pattern=custom north west lines,
            hatchspread=20pt, hatchthickness=7pt, hatchshift=5pt,
            minimum width=\textwidth + 0.375in, minimum height=1.75in,
            anchor=center]
            (linesrect) at ($(goldrect.center) + (0in,0in)$) {};
      \node[anchor=center, text width=\textwidth, outer sep=0in, inner sep=0.125in,
            text shadow={[align=center] at (goldrect.center) {\Huge \champion #1\\ \LARGE \champion #2}}]
      at (goldrect.center)
            {\Huge \champion \color{grey60} #1\\ \LARGE \champion #2};
    \end{pgfonlayer}
  \end{tikzpicture}%
}

\newlength{\cellheight}
\newlength{\cellwidth}
\NewEnviron{pucell}[4]
{
  \noindent%
  \ifnum\incell=1
    \setlength{\cellheight}{1.875in * \real{#2} - 0.1in * \real{#2} * \real{\numheadings}}%
    \setlength{\cellwidth}{3.625in * \real{#1}}%
  \else
    \setlength{\cellheight}{2.0in * \real{#2} - 0.09in * \real{#2} * \real{\numheadings}}%
    \setlength{\cellwidth}{3.75in * \real{#1}}%
  \fi
  %\begin{tikzpicture}[fill opacity=.4,remember picture, overlay]
  %  \fill[red!40] (0,0) rectangle (\cellwidth, -\cellheight);
  %\end{tikzpicture}%
  \begin{minipage}[t][\cellheight][t]{\cellwidth}
    \noindent%
    % Check if we're already inside a cell, if not, we make the outer title
    \ifnum\incell=1
      \ifx&#3&%
        % #1 is empty
      \else
        %\hspace*{-0.625in}%
        \innertitle{#3}{#4}%
        \vspace{1.75in}
      \fi
    \else
      \outertitle{#3}%
      \vspace{1.25in}\\
    \fi
    \def\incell{1}%
    %\vspace*{\fill}
    \BODY
    %\vspace*{\fill}
  \end{minipage}%
}

\NewEnviron{pucol}[2]
{
  \setlength{\cellheight}{1.875in * \real{#2} - 0.083in * \real{#2} * \real{\numheadings}}%
  \setlength{\cellwidth}{3.625in * \real{#1}}%
  \begin{minipage}[t][\cellheight][t]{\cellwidth}
    \BODY
  \end{minipage}
}

\newcommand{\posterfigure}[3]{
  \filename@parse{#1}
  \begin{figure}[H]
    \ifthenelse{\equal{\filename@ext}{pgf} \OR \equal{\filename@ext}{tex}}{
      \centering\adjustbox{keepaspectratio=true, height=\cellheight-4em}{\input{#1}}
    }{
      \centering\includegraphics[keepaspectratio=true, width=\cellwidth, height=\cellheight]{#1}
    }
  \centering \caption{#2}\label{#3}
  \end{figure}%
}

\def\labelitemi{\Huge $\Rightarrow$}
