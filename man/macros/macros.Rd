\newcommand{\CO2}{\ifelse{latex}{\eqn{\mathrm{CO_{2}}}}{\ifelse{html}{\out{CO<sub>2</sub>}}{CO2}}}
\newcommand{\NH3}{\ifelse{latex}{\eqn{\mathrm{NH_{3}}}}{\ifelse{html}{\out{NH<sub>3</sub>}}{NH3}}}
\newcommand{\H2S}{\ifelse{latex}{\eqn{\mathrm{H_{2}S}}}{\ifelse{html}{\out{H<sub>2</sub>S}}{H2S}}}
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}
\newcommand{\O2}{\ifelse{latex}{\eqn{\mathrm{O_{2}}}}{\ifelse{html}{\out{O<sub>2</sub>}}{O2}}}

\newcommand{\ZC}{\ifelse{latex}{\eqn{Z_\mathrm{C}}}{\ifelse{html}{\out{<I>Z</I><sub>C</sub>}}{ZC}}}
% use \nH2O{&#772;} to call this macro (the html code can't be defined in the macro,
% which interprets '#' followed by a number as a placeholder for an argument)
\newcommand{\nH2O}{\ifelse{latex}{\eqn{\bar{n}_{\mathrm{H_{2}O}}}}{\ifelse{html}{\out{<i>n</i>#1<sub>H<sub>2</sub>O</sub>}}{nH2O}}}

\newcommand{\logfO2}{\ifelse{latex}{\eqn{\log f_{\mathrm{O_{2}}}}}{\ifelse{html}{\out{log<i>f</i><sub>O<sub>2</sub></sub>}}{logfO2}}}
\newcommand{\logaH2O}{\ifelse{latex}{\eqn{\log a_{\mathrm{H_{2}O}}}}{\ifelse{html}{\out{log<i>a</i><sub>H<sub>2</sub>O</sub>}}{logaH2O}}}
