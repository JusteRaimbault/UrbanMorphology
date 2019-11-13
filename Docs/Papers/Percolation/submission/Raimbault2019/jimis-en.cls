\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{jimis-en}

\LoadClass[a4paper,12pt]{article}

\RequirePackage[a4paper,left=2.5cm,right=2.5cm,top=2.5cm,bottom=2.7cm]{geometry}
\RequirePackage{times}
\RequirePackage{authblk}
\renewcommand{\Authands}{, }
\renewcommand{\Authfont}{\bfseries}

\RequirePackage[usenames,dvipsnames,table]{xcolor}
%\definecolor{epicolor}{RGB}{14,159,149}
\definecolor{epicolor}{RGB}{12,135,127}

\RequirePackage[
    bookmarksnumbered=true,bookmarksopen=true,
	unicode=true,colorlinks=true,linktoc=all,%linktoc=page
	linkcolor=epicolor,citecolor=epicolor,filecolor=epicolor,urlcolor=epicolor,
	pdfstartview=FitH
]{hyperref}

\RequirePackage{parskip}
\RequirePackage{graphicx}
\RequirePackage{fancyhdr}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FROM A PREPRINT TO A PAPER
%% PPREPRINT
% \lfoot{\parbox[t]{0.5\textwidth}{\scriptsize  
%   {} 
%   }}
% \rfoot{\renewcommand{\UrlFont}{\normalfont\bfseries\scriptsize} 
%   {} 
% }
% \rhead{}
% \lhead{}
% \renewcommand{\headrulewidth}{0pt}
% \pagestyle{fancy}
%% PAPER
\lfoot{\parbox[t]{0.5\textwidth}{\scriptsize  
  \textbf{J. of Interd. Method. and Issues in Science} \\
  Open-access journal: \url{http://jimis.episciences.org}
  }}
\rfoot{\renewcommand{\UrlFont}{\normalfont\bfseries\scriptsize} 
  \scriptsize{%ISSN: xxxxxxx, 
  \copyright JIMIS, Creative Commons} \\
  \scriptsize{\@publication, \@doi} 
}
\rhead{}
\lhead{}
\renewcommand{\headrulewidth}{0pt}
\pagestyle{fancy}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[largestsep]{titlesec}
\titleformat{\section}{\normalfont\normalsize\bfseries}{\thesection}{1em}{\MakeUppercase}
\titleformat{\subsection}{\normalfont\normalsize\bfseries}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\normalfont\normalsize\itshape}{\thesubsubsection}{1em}{}
\titleformat{\paragraph}[runin]{\normalfont\normalsize\bfseries}{\theparagraph}{1em}{}[. \mbox{}]
\titleformat{\subparagraph}[runin]{\normalfont\normalsize\itshape}{\thesubparagraph}{1em}{}[. \mbox{}]

\titlespacing*{\section}      {0pt}{1.5ex plus .7ex minus .2ex}{-\parskip \advance\parskip by 0pt plus .2ex}
\titlespacing*{\subsection}   {0pt}{1.5ex plus .7ex minus .2ex}{-\parskip \advance\parskip by 0pt plus .2ex}
\titlespacing*{\subsubsection}{0pt}{1.5ex plus .7ex minus .2ex}{-\parskip}
\titlespacing*{\paragraph}    {0pt}{.8ex plus .5ex minus .2ex}{1em}
\titlespacing*{\subparagraph} {0pt}{.8ex plus .5ex minus .2ex}{1em}

\renewcommand{\thesection}{\Roman{section}}
\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}

\RequirePackage{float} % must be loaded before caption package
\floatstyle{ruled}
\newfloat{listing}{tbph}{lol}[section]
\floatname{listing}{Algorithme}
\renewcommand{\floatpagefraction}{0.8}
\renewcommand{\bottomfraction}{0.5}
\renewcommand{\topfraction}{0.5}
\renewcommand{\textfraction}{0.1}

\RequirePackage{rotating,subfigure} % must be loaded before caption package
\RequirePackage[small]{caption}
\setlength{\abovecaptionskip}{1ex plus 0.5ex minus 0.2ex}

\RequirePackage[fleqn]{amsmath}
\RequirePackage{amsthm}
\newtheoremstyle{jimis}% name of the style to be used
  {2ex plus 0.2ex minus 0.2ex}% measure of space to leave above the theorem. E.g.: 3pt
  {1ex plus 0.2ex minus 0.2ex}% measure of space to leave below the theorem. E.g.: 3pt
  {\normalfont}% name of font to use in the body of the theorem
  {0pt}% measure of space to indent
  {\itshape}% name of head font
  {\\}% punctuation between head and body
  {0pt}% space after theorem head; " " = normal interword space
  {\thmname{#1}\thmnumber{ #2}: \thmnote{ #3}}% Manually specify head
\theoremstyle{jimis}
\newtheorem{definition}{D\'efinition}

\RequirePackage{enumitem}
\setlist{noitemsep,topsep=-1ex}

\RequirePackage{listings}
\lstset{numbers=left,numberstyle=\scriptsize,basicstyle=\small,xleftmargin=3ex,flexiblecolumns=true}

\RequirePackage{natbib}
\renewcommand\bibsection{\subsection*{\refname}}
\renewcommand{\bibfont}{\footnotesize}
\setlength{\bibsep}{5pt}

% ------------------------------------------------------------- %

% corresponding author
\newcommand{\@corrauthor}{}
\newcommand{\corrauthor}[2]{\renewcommand{\@corrauthor}{*Corresponding author: \href{mailto:}{#1} \if\relax#2\relax\else (\texttt{#2}) \fi}}

%DOI number
\newcommand{\@doi}{}
\newcommand{\doi}[2]{\renewcommand{\@doi}
{\textit{DOI: \href{https://doi.org/#1}{#1}} \if\relax#2\relax\else- \textit{#2} \fi}}

% Review process with submission and acceptance dates
\newcommand{\@review}{}
\newcommand{\review}[2]{\renewcommand{\@review}
{\textit{Submitted: {#1}} \if\relax#2\relax\else- \textit{Published: {#2}} \fi}}

% Volume and Year
\newcommand{\@publication}{}
\newcommand{\publication}[2]{\renewcommand{\@publication}
{Volume: \textbf{#1} \if\relax#2\relax\else- Year: \textbf{#2}\fi}}

% Issue 
\newcommand{\@issue}{}
\newcommand{\issue}[1]{\renewcommand{\@issue}
{Issue: \textbf{#1}}}

% Editors
\newcommand{\@editors}{}
\newcommand{\editors}[1]{\renewcommand{\@editors}
{Editor: \textit{#1}}}

% maketitle identical to standard article title, except for pagestyle
\renewcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    %\thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}

\renewcommand{\@maketitle}{%
  \newpage
  \null
  \vspace{-2cm}
  \includegraphics[height=1.7cm]{episcience_logo.png}\hfill 
  \includegraphics[height=3cm]{jimis_logo.png}\hfill 
  \includegraphics[height=1.6cm]{espace_logo.png}
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\large\bfseries \@title \par}%
    \vskip 1.5em%
    {\small
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author\\[\affilsep]
        \@corrauthor \\ \\
        \@doi \\
        \@review \\
        \@publication \\
        \@issue \\
        \@editors
%        \vspace{0.3cm}
      \end{tabular}\par}%
  \end{center}%
  \par
  \vskip 1.5em}

\renewcommand{\abstract}[1]{%
  \hrule
  \vspace{9pt}%
  \begingroup
  \leftskip0.5cm
  \rightskip\leftskip
  {\small
  {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}\par%
  #1}%
  \par
  \endgroup}

\newcommand{\keywords}[1]{%
  \begingroup
  \leftskip0.5cm
  \rightskip\leftskip
  \par\vspace{1ex}{\small
  {\bfseries Keywords\vspace{-.5em}\vspace{\z@}}\par%
% {\bfseries Mots-Cl\'es \vspace{-.5em}\vspace{\z@}}\par%
  #1}%
  \vspace{14pt}%
  \hrule
  \vspace{10pt}%
  \par
  \endgroup}