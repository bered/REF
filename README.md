
        <h1 align=center>REF - Rational Ess Finder</h1>
        <h2 align=center>A solver for Standard Quadratic Problems</h2>

        <h2>Introduction</h2>
        <p>REF is a command-line C++ program for calculating all strict local maximizers of a Standard Quadratic Problem (StQP). From the point of view of Evolutionary Game Theory this is equivalent to finding all Evolutionary Stable Strategies (ESSs) of a partnership game. The algorithm for doing so was originally described in I.M. Bomze: "Detecting all evolutionarily stable strategies", Journal of Optimization Theory and Applications 75.2 (1992), and this implementation is a further development and  improvement of it. It was most recently used <a href="http://www.optimization-online.org/DB_HTML/2016/05/5452.html" target="_blank">here</a>, although the implementational details are not explained.</p>

        <p>REF takes a rational, not necessarily symmetric n &times n-matrix as an argument and outputs the number of ESSs (as well as detailed information about them) in exact arithmetics. It does so by traversing (potentially) all the 2^n support sets of {1,...,n}, and REF can (theoretically) handle matrices up to n=64. In practice, depending on your CPU, your memory and your time, n up to 30 or 40 is feasible.</p>

        <h2>Downloads</h2>
        <ul>
            <li><a href="REF_win64.zip" target="_blank">Windows 64bit</a></li>
            <li><a href="REF_ubuntu64.zip" target="_blank">Ubuntu/Debian 64bit</a></li>
            <li>MacOsX - cross-compilation is quite complicated, please compile it yourself!</li>
            <li><a href="REF_source.zip" target="_blank">C++ Source</a> (as a Code::Blocks Project)</li>
        </ul>

        <h2>Initial Test (for Windows)</h2>
        <ul>
            <li>Download the .zip under "Windows 64bit".</li>
            <li>Create a folder "C:\REF" and unpack the .zip there.</li>
            <li>Initial test: open a command line and copy the following code: <samp>C:\REF\REF.exe 21#15,15,7,15,15,7,7,15,7,15</samp></li>
            <li>REF should run for 5-30 seconds (depending on your hardware) and write "4410" as output.</li>
        </ul>

        <h2>Documentation</h2>
        Now you are able to use REF, this is possible in two ways. <br>
        <ul>
            <li>You manually type into the command line - as done for the initial test. This is only recommended for experimenting with the command line options or for quick checks.</li>
            <li>You call the command line from a script/program of your choice. This is the way the program is designed to run.</li>
        </ul>
        The program can be seperated in three parts: <samp>REF.exe [options] [matrix]</samp> <br>
        <ul>
            <li><samp>REF.exe</samp> is just the name of the program you call, it is always the same.</li>
            <li><samp>[options]</samp> are the so-called command-line-options, see below.</li>
            <li><samp>[matrix]</samp> is the matrix you want to search for ESS, see below.</li>      
        </ul>

        <h3>Input [options]</h3>
        <ul>
            <li><strong>-v or --vectors</strong> <br>
                Usually the program only outputs the number of ESS a matrix has, this is useful if you want to search many matrices fast. The option <samp>-v</samp> enables you to output the whole information which is stored about the solution of your problem. This is done in a .csv-style manner, see section "Resulting Vectors".</li>
            <li><strong>-f or --fullsupport</strong> <br>
                Usually the program starts with the smalles support sets and searches up to the bigger ones. This behaviour is favoured if you expect the matrix to have many ESS. But if you expect the matrix to have only one ESS with full support, you should set this option. Now first the support sets of size one are searched, next the whole support is searched, and afterwards the support sets from size 2 to n-1 are searched in increasing order.</li>
            <li><strong>-e or --exact</strong> <br>
                Usually most of the calculations are done as floating points (double), and only if the outcome seems to be a solution then the calculation is repeated with rational numbers. This vastly increases the speed of the algorithm and works in almost all cases. But if the differences of the magnitude of the matrix inputs are huge (say more than 10^10), this could lead to unnoticed errors, since the double-datatyp only has a precision of 15-17 digits. This option disables the use of floating points, all calculations are done entirely as rationals. You can also enable this option to double-check a result.</li>
            <li><strong>-l or --log</strong> <br>
                This option generates a detailed log-file named "Log.txt" in the directory where the program runs. For diagnostic- or learning-purposes only!</li>
        </ul>

        <h3>Input [matrix]</h3>
        At the moment two different ways of passing a matrix are supported
        <ul>
            <li><strong>Whole matrix</strong> <br>
                Example: <samp>2#5,6/7,-1,-2/3</samp> <br>
                First the dimension of the matrix is given, then a hash-symbol (#), afterwards the numbers as integers or rationals, where rationals are written with a slash (/). All numbers can be positive (nothing) or negative (-) and are seperated by a comma (,). Whitespaces or other characters are not allowed! The numbers must be in row-major-order, that means that the first <em>row</em> of the matrix is given as (5,6/7) and the second row as (-1,-2/3).</li>
            <li><strong>Cyclically symmetric matrix</strong> <br>
                Example: <samp>7#2/3,5,9</samp> <br>
                Since cyclically symmetric matrices play a major role in the search for hard cases for StQP's, there is an shortcut for inputting them. The input in the example generates a 7x7 matrix where the first line is given as (2/3,5,9,9,5,2/3,0), ie. take the vecor, glue it together with its mirrored version and glue a zero at the end. Now this row is the last row of the matrix, all other rows are generated by rotational left shifts of order n-i, where i is current row index. Note that if n is even then the last number of the input is not doubled. Use this type of input also if you want to enable the optimziations for cyclically symmetric matrices, see section "Cyclically Symmetric Opimization".</li>
        </ul>

        <h3>Full Examples</h3>
        <samp>REF.exe -v 2#5,6/7,-1,-2/3</samp> <br>
        <samp>REF.exe -l -f -e 19#1,2,2,2,2,2,1,1,2</samp> <br>
        <samp>REF.exe 23#27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939</samp>
        <p>A script in R which tests some matrices in multi-threading can be found <a href="ref_tester.R" target="_blank">here</a>.</p>

        <h3>Resulting Vectors</h3>
        If the option -v is enabled, then the lines of output of REF consist of the following: <br>
        <samp>[Number of ESSs]<br>[CSV Header]<br>[CSV Data]<br>[CSV Data]<br>[CSV Data]<br>...</samp>
        <p>One easy way to process the data is to seperate the first from the following lines, and parse these lines by a .csv-parser. The field sperator is a semicolon, field delimiters are not used.</p>
        See a full example here, the output of <samp>~/REF -v 5#1,0,2,2,2,0,1,2,2,2,2,2,1,0,0,2,2,0,1,0,2,2,0,0,0</samp> (input under Ubuntu): <br>
        <samp>
            4<br>
            VectorID;Vector;Support;SupportSize;ExtendedSupport;ExtendedSupportSize;ShiftReference;IsEss;Reason;Payoff;PayoffDecimal<br>
            1;1/2,0,1/2,0,0;5;2;5;2;0;1;2;3/2;1.500000<br>
            2;0,1/2,1/2,0,0;6;2;6;2;0;1;2;3/2;1.500000<br>
            3;1/2,0,0,1/2,0;9;2;9;2;0;1;2;3/2;1.500000<br>
            4;0,1/2,0,1/2,0;10;2;10;2;0;1;2;3/2;1.500000<br>
            5;2/3,0,0,0,1/3;17;2;29;4;0;0;7;4/3;1.333333<br>
            6;0,2/3,0,0,1/3;18;2;30;4;0;0;7;4/3;1.333333<br>
        </samp>
        <p>REF with option -v does not only save all the ESSs, but all the "Candidates". A serious candidate for an ESS is an isolated Nash Equilibrium strategy (equivalently an isolated KKT-point), which is in the relative interior of the considered support, but only if no other candidate's support is contained in its support. Sloppily speaking, only the smaller isolated Nash Equilibrium strategies are candidates. For details see the references in the "Introduction".</p>
        Description of the different fields:<br>
        <ul>
            <li><strong>VectorID</strong> <br>
                Auto-incrementing integer.</li>
            <li><strong>Vector</strong> <br>
                Returns the vector of the ESS of dimension n as a rational number. The encoding is the same as for the input matrices.</li>
            <li><strong>Support</strong> <br>
                The support set of the ESS coded in binary, e.g. support 6 means that the support set is {2,3}. Always smaller than 2^n.</li>
            <li><strong>SupportSize</strong> <br>
                The amount of numbers in the support set</li>
            <li><strong>ExtendedSupport</strong> <br>
                The extended support coded in binary. See the references in the "Introduction" for a definition.</li>
            <li><strong>ExtendedSupportSize</strong> <br>
                The amount of numbers in the extended support set.</li>
            <li><strong>ShiftReference</strong> <br>
                Only used for cyclically symmetric matrices, see section "Cyclically Symmetric Optimization", otherwise 0.</li>
            <li><strong>IsEss</strong> <br>
                Is 1 if the candidate is an ESS, otherwise 0.</li>
            <li><strong>Reason</strong> <br>
                Since the procdure to determine if a candidate is an ESS is quite cumbersome, here the reason for the decision is stated as an integer. The coding is the following: <br>
                true_pure_ess=1<br>
                true_posdef_double = 2<br>
                true_posdef_rational=3<br>
                true_copositive=4<br>
                false_not_posdef_and_K_0_1 = 5<br>
                false_not_partial_copositive = 6<br>
                false_not_copositive=7<br>
                Numbers 1 to 4 mean that the candidate is an ESS, 5 to 7 mean that it is not. For details see the references in the "Introduction" and the source code.</li>
            <li><strong>Payoff</strong> <br>
                The value of the maximum for this particular maximizer (or equivalently the average population payoff for the ESS) as a rational number.</li>
            <li><strong>PayoffDecimal</strong> <br>
                The same as a decimal number.</li>
        </ul>

        <h3>Cyclically Symmetric Opimization</h3>
        Whenever the shortcut for inputting cyclically symmetric matrices is used, a certain optimization is performed, utilizing the symmetrie properties of the matrix. Careful, if a cyclically symmetric matrix is inputted without the shortcut, then this optimization is not performed! <br>
        If a candidate/ESS for a cyclically symmetric matrix is found and the support size of the candidate/ESS and n are coprime, then it is known that another n-1 candidates/ESSs exist, they are just shifted (rotated) versions of the original candidate/ESS. The optimization exploits that fact, and the n-1 other candidates/ESSs are just recorded and not calculated. This increases the speed of REF, depending on n and the size of the numbers of the input matrix. For n being prime the increase is clearly the most. <br>
        The field ShiftReference for the vector-output references the VectorID of the originally calculated ESS (and this ESS also references itself). For analyzing this cyclical structure this can be used by querying the .csv for VectorID==ShiftReference.

        <h2>License</h2>
        The software REF is open-source under the GNU General Public License, Version 3. No warranty or liability of any kind is provided.<br>
        Other open-source software used by REF, with special thanks:
        <ul>
            <li><a href="http://www.boost.org/" target="_blank">Boost C++ libraries</a></li>
            <li><a href="http://www.argtable.org/" target="_blank">Argtable3</a></li>
        </ul>

        <h2>Contact</h2>
        For questions, bug reports and suggestions please write to: reinhard.ullrich (at) univie.ac.at.<br>
        Last updated 19.11.2016.
