\NeedsTeXFormat{LaTeX2e}[1998/06/01]
\ProvidesClass{puposter}[2016/06/22 class for purdue poster design]

%% Load base
\LoadClass{article}

\RequirePackage{fontspec}
\RequirePackage{tikz} % to draw
\RequirePackage{anyfontsize} % to use the fonts
\RequirePackage{lipsum}
\RequirePackage{etoolbox}
\RequirePackage{pst-all}
\RequirePackage{pst-light3d}
\RequirePackage{lmodern}
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

% defining the new dimensions and parameters
\newlength{\hatchspread}
\newlength{\hatchthickness}
\newlength{\hatchshift}
\newcommand{\hatchcolor}{}
% declaring the keys in tikz
\tikzset{hatchspread/.code={\setlength{\hatchspread}{#1}},
         hatchthickness/.code={\setlength{\hatchthickness}{#1}},
         hatchshift/.code={\setlength{\hatchshift}{#1}},% must be >= 0
         hatchcolor/.code={\renewcommand{\hatchcolor}{#1}}}
% setting the default values
\tikzset{hatchspread=3pt,
         hatchthickness=0.4pt,
         hatchshift=0pt,% must be >= 0
         hatchcolor=black}
% declaring the pattern
\pgfdeclarepatternformonly[\hatchspread,\hatchthickness,\hatchshift,\hatchcolor]% variables
   {custom north west lines}% name
   {\pgfqpoint{\dimexpr-2\hatchthickness}{\dimexpr-2\hatchthickness}}% lower left corner
   {\pgfqpoint{\dimexpr\hatchspread+2\hatchthickness}{\dimexpr\hatchspread+2\hatchthickness}}% upper right corner
   {\pgfqpoint{\dimexpr\hatchspread}{\dimexpr\hatchspread}}% tile size
   {% shape description
    \pgfsetlinewidth{\hatchthickness}
    \pgfpathmoveto{\pgfqpoint{0pt}{\dimexpr\hatchspread+\hatchshift}}
    \pgfpathlineto{\pgfqpoint{\dimexpr\hatchspread+0.15pt+\hatchshift}{-0.15pt}}
    \ifdim \hatchshift > 0pt
      \pgfpathmoveto{\pgfqpoint{0pt}{\hatchshift}}
      \pgfpathlineto{\pgfqpoint{\dimexpr0.15pt+\hatchshift}{-0.15pt}}
    \fi
    \pgfsetstrokecolor{\hatchcolor}
%    \pgfsetdash{{1pt}{1pt}}{0pt}% dashing cannot work correctly in all situation this way
    \pgfusepath{stroke}
   }
\pgfdeclarepatternformonly[\hatchspread,\hatchthickness,\hatchshift,\hatchcolor]% variables
  {custom north east lines}% name
  {\pgfqpoint{\dimexpr-2\hatchthickness}{\dimexpr-2\hatchthickness}}% lower left corner
  {\pgfqpoint{\dimexpr\hatchspread+2\hatchthickness}{\dimexpr\hatchspread+2\hatchthickness}}% upper right corner
  {\pgfqpoint{\dimexpr\hatchspread}{\dimexpr\hatchspread}}% tile size
  {% shape description
   \pgfsetlinewidth{\hatchthickness}
   \pgfpathmoveto{\pgfqpoint{\dimexpr\hatchshift-0.15pt}{-0.15pt}}
   \pgfpathlineto{\pgfqpoint{\dimexpr\hatchspread+0.15pt}{\dimexpr\hatchspread-\hatchshift+0.15pt}}
   \ifdim \hatchshift > 0pt
     \pgfpathmoveto{\pgfqpoint{-0.15pt}{\dimexpr\hatchspread-\hatchshift-0.15pt}}
     \pgfpathlineto{\pgfqpoint{\dimexpr\hatchshift+0.15pt}{\dimexpr\hatchspread+0.15pt}}
   \fi
   \pgfsetstrokecolor{\hatchcolor}
%    \pgfsetdash{{1pt}{1pt}}{0pt}% dashing cannot work correctly in all situation this way
   \pgfusepath{stroke}
  }

\pagestyle{empty}

\renewcommand\maketitle{
    \begin{tikzpicture}[remember picture, overlay]
        \begin{pgfonlayer}{lines}
            \fill[newgold]
                  ($(current page.north west) + (1.25in, -1.25in)$)
                  rectangle
                  ($(current page.north east) + (-0.75in, -4.25in)$);
            \fill[hatchcolor=white, pattern=custom north west lines,
                  hatchspread=30pt, hatchthickness=10pt]
                  ($(current page.north west) + (1.25in, -1.25in)$)
                  rectangle
                  ($(current page.north east) + (-0.75in, -4.50in)$);

            \fill[grey60]
                ($(current page.north west) + (1in, -1in)$) --
                ($(current page.north east) + (-1in, -1in)$) --
                ($(current page.north east) + (-1in, -4in)$) --
                ($(current page.north west) + (1in, -4in)$) --
                cycle;
        \end{pgfonlayer}
        \begin{pgfonlayer}{main}
            \node[
              align=left,
              font=\champion\color{newgold}\fontsize{72}{80}\bfseries,
              anchor=north west,
              minimum width=\paperwidth - 2.75in,
              text width=\paperwidth - 2.75in, inner sep=0, outer sep=0
              ]
              at ($(current page.north west) + (1.375in, -1.375in)$)
              (title)
              {\@title};
        \end{pgfonlayer}
    \end{tikzpicture}
}

\AtBeginDocument{
\begin{tikzpicture}[remember picture,overlay]
\begin{pgfonlayer}{background}
\node[scale=1.5, inner sep=0, outer sep=0,
      anchor=south east]
     (bottomright) at (current page.south east)
     {\includegraphics{css/background_bottom_right.pdf}};
\node[scale=1.5, inner sep=0, outer sep=0, anchor=north west] (topleft) at
    (current page.north) {\includegraphics{css/background_top_left.pdf}};
\end{pgfonlayer}

\end{tikzpicture}
}
