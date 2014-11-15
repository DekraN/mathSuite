mathSuite
//-----------------------------------------------------------------------------------------------------------------//
AUTHOR: Marco Chiarelli aka DekraN aka Wesker013 (FB)								   ||
CURRENT VERSION: v7.00											   	   ||
LAST UPDATE: 16:00 15/11/2014											   ||
CONTACT ME at: marco_chiarelli@yahoo.it										   ||
or at marcochiarelli.nextgenlab@gmail.com									   ||
//-----------------------------------------------------------------------------------------------------------------//
///CHANGELOG v2.5: Fixed every imaginable bug and started     							  ///
///    object-oriented style programming in the code skeleton, by treating the programs structured variables      ///
///                       as an object-oriented-language class definitions.                                       ///
/// Introduced Families Programs System, as a new type, which contains other childs programs, and one of these    ///
/// could be itself a new access point for a new Programs Family. Sincerely don't know whether this system could  ///
///     improve or not suite performances. I'll ask more informations about this to the Professor Bochicchio.     ///
///     Added new simplifying Vectorial MACRO-Operations and Introduced Inverse Operations for L.A. Programs.     ///
///             You could use them by enabling associated setting in the Settings Program.                        ///
/// If the Automatic Deactivator for the I.O. is Active, then the program turns back to Normal Operations after   ///
/// each Opposite Operation Process. Else, Opposite Operations will remain active until you disable them from SM  ///
///                                      CHANGELOG v3.0 (30/04/2013) -                                            ///
/// Added rank calculator for a [nXm] Matrix, inverse calculator for a Square Matrix, transpose calculator for a  ///
/// [nXm] Matrix. Now it's possible to set some default values and pre-compiling params from CONFIG.INI, whose    ///
/// environment I've made the most user-friendly possible. Completely redefined the code and increased to more or ///
/// less 15% calculus speed level of various subprograms. Enhanced explicit cast performance in calcolatoreDiBase ///
///   subprogram, and Added char_insert boolean value in order to check whether user inserts char or not in the   ///
///                      various Matrix inserting checkpoints of the entire program.                              ///
///                                      CHANGELOG v3.2 (09/05/2013) -                                            ///
/// Added Dynamic Memory Allocation concerning Matrix. Automatic Matrix Inserting System is now more simple and   ///
/// faster. Allocation System is M.L.R.F. (Memory Leak Risks Free), because I checked every single instruction    ///
/// concerning enterMatrix(...) function and resulted 100% protected. There is no risk malloc, realloc and free   ///
///   functions used by me could enter in another memzone except the one and only my program needs. Introduced    ///
///            (Initializer+Stabilizer) for Matrix Allocation. Normally my program would allocate every           ///
/// new columns or raws's necessary size every time a new matrix item is entered by A.M.I.S. (upwriten meaning)   ///
/// If Initializer+Stabilizer binary setting is enabled, then INITIALIZER will allocate a quantum of memory size  ///
/// necessary for various both raws and columns BEFORE that User could start to Enter The Matrix. This memory     ///
/// size depends by 'ityp' program default type size, and by a constant value first defined in CONFIG.INI         ///
/// but changeable during program execution), named fattore_stabilizzatore (Default Value) or STABILIZER_FACTOR,  ///
/// which is the variable name that is physically located into memory. Every time User trespasses this memorysiz  ///
/// quantum, INITIALIZER+STABILIZER action will provide to autoincrease its size, exactly like the way            ///
/// updescripted. At the end of the User INPUT, STABILIZER will reallocate the matrix, by eliminating the unused  ///
///   memory which INITIALIZER+STABILIZER allocated for the last time. Don't know whether this memory is faster   ///
///      than Normal Mode of Dynamic Memory Matrix Allocation or not. As always, I'll discuss this problem with   ///
/// Professor Bochicchio ASAP. My first hypothesis is that Normal Mode is more performing, because S+I method     ///
/// will give less instructions to the CPU to process, but HEAP memory region size will lineary increase          ///
/// depending by STABILIZER_FACTOR; and so I suppose you know, memory is not faster as CPU. Now Matrix Inserting  ///
/// system could be defined in CONFIG.INI and also changed during program execution. The two possibilities are    ///
/// Automatic (normally used) and non-Authomatic, that requests to the User the Matrix Dimensions every time he   ///
///     have to Enter a Matrix. Now it's possible to set Program Default Type, by modifying the relative ITAnamed ///
/// variable onto CONFIG.INI. This type could be float, double, long double, int, short, long int, long long int  ///
///   and also char!!! This means you could work with Chars Matrix, and it's not a useless feature, since C char  ///
///                                  type is native lexicographic ordered.                                        ///
///                                      CHANGELOG v3.35 (11/05/2013) -                                           ///
///      Introduced famous successions calculus, like Fibonacci, Geometric and the Generalized Armonic ones       ///
///                  and also the possibility of calculate their sum at a specified INDEX.                        ///
///      Introduced Arythmetic and Geometric Media calculus, and Added simple Summation System.                   ///
///          All these features are simply achievable by starting calcolatoreDiBase subprogram                    ///
///                                      CHANGELOG v3.5 (15/05/2013)                                              ///
/// Removed all type strict controls. Now, even if ityp default type is pre-compilation changeable, whereas it's  ///
/// directly involved in INPUT/OUTPUT FORMATTATIONS, it will be necessary some modifies in this script zones each ///
/// time you change it. However, I put in CONFIG.INI some facilitations about this not simply work. Introduced    ///
/// OVERFLOW CHECKER in all INPUT/OUTPUT numeric elaborations. If Domain Check mode is enabled, then System will  ///
/// print you the classic "OVERFLOW ERROR". Fixed a bug into DivisionModes-0-Inserting-Proof System. Now it works ///
/// perfectly. Code cleaned and highly improved executing speed. I also inlined some basic not-complex functions. ///
/// Obviously this had the bad inevitable side-effect of increasing total script size. But it's a good affair     ///
/// since nowaday progress difference between CPUs and memory is unimaginable.                                    ///
/// Improved media and media_geometrica, which are the two relevant features added in the last version in Base    ///
///     Calculator subprogram. Now both requires the b parameter to be entered. If it isn't 0, then Media will    ///
///   increase, else it don't. Improved all the inlined prestigious successions functions like fibonacci and      ///
///       factorial, by adding unsigned INVERSE_OPS at the int64_t (long long int) formal parameters type.        ///
///              Introduced Armonic Mean, Power Mean, Mean Value and Median in Base Calculator subprogram.        ///
///                                      CHANGELOG v4.0 (23/05/2013)                                              ///
/// Introduced new config.INI Initial Settings Configurations, but by excluding Back to MAIN MENU Char and all    ///
/// the Operations Chars, because GetPrivateProfileString is giving me problem, both Unicode (W), both ANSI (A)   ///
/// version, tested respectively on UNICODE config.INI and ANSI config.INI. However, I let Standard ANSI config   ///
/// INI file in the Official Version. Moved all DEFAULTS Program Settings Value Constants into defaults.h new     ///
/// header, correctly linked to the Current Project. Introduced new Showing (Sub)Programs Descriptions by reading ///
/// the correspondent file, whose name is formatted following a certain algorythm rule that is defined into       ///
/// engageDescriptions function, located in geometry.c library file. Introduced MsgBox system for Long Strings    ///
/// You could try its usefulness by requesting to the program to Show a certain Program Description. This was     ///
/// implemented for not to fill excessively the stdout buffer, whose dimensions are very important for a fast and ///
/// clean execution of the calculus environment. Transformed the Program Information STATIC (SUB)MASK, whose name ///
/// in the program is "Informazioni sul PROGRAMMA", into a Dynamic Macro, where now the relative Informations     ///
/// it gives dinamically refers to the Current Father (Sub)Program Description. Fixed a bug in the matrixProduct  ///
/// function, which was causing an HEAP OVERFLOW with the Lazy Execution bool value disabled (now I've renamed    ///
/// STABILIZER+INITIALIZER System's correspondent bool value name in Lazy Execution. If the Lazy Execution is     ///
/// enabled, then STABILIZER+INITIALIZER action is not present, else there will be the mechanism here descripted  ///
/// concerning Dynamic Matrixes Mallocation. This is an evident tribute to Lazy Evaluation Programming            ///
/// Languages technique of evaluating a statement or an expressions possibily only when It is requested           ///
/// effectively. Performed the ERRORS Management System by using external variables errno, located into errno.h   ///
///native header file. Now you could gain more informations on happening of a certain Error, by viewing the ERROR ///
///                                        TAG/ID that will be printed.                                           ///
///                                      CHANGELOG v4.2 (26/05/2013)                                              ///
/// Optimized Code and Fixed some bugs into matrixProduct, precisely in the MATRIX_POWER subProgram. Now the      ///
/// product system is more efficient and readable. Now MATRIX_POWER gives in OUTPUT correct results.              ///
/// Enhanced Domain Check System into Trigonometric Operations, by inserting a new inline function which checks   ///
/// by using a loop whether the number is equal to M_PI_2 + k*M_PI. This may be useful for example to check if    ///
/// a certain INPUT Parameter belongs to the tan(x) Domain. The Precision of this System is not guaranteed if     ///
/// use Radiant Parameters Entering System. Yeah, you could now choice your Trigonometric Operations Parameters   ///
/// Entering Modality, by managing its own boolean settings, located into Change PROGRAM Settings SubProgram.     ///
/// If you want to use a High-Level TRIGONOMETRIC_DOMAINs Checker, you have to enable Degrees Parameters Entering ///
/// System, because the Degrees to Radiant Conversion System, gives the System a more precise Result to the       ///
/// various functions, and in this case, to the trigonometric_domain(n) inline function, both declared and        ///
/// defined in a new Header, called ext_math.h as EXTERNAL Math Library, where some constants and functions has   ///
/// been re-defined. In this Project Header I put all the Properly Math Functions, both declarations and both     ///
///                              definitions, for increasing the Code Readability.                                ///
///          cot(n), csch(n), sech(n), coth(n), acsc(n), asec(n), acot(n), acsch(n), asech(n), acoth(n);          ///
///   For More Informations, see MAT-LAB Calculus Environment infos, in which you can find the Full Documentation ///
///                                      about their definitions, domains, etc.                                   ///
/// Added also cbrt(...) Standard Function into Base Calculator subProgram, which calculates the Cubic Root of a  ///
/// given ityp INPUT Parameter. Enhanced the calcolatoreDiBase() related function, particularly the Operations ID ///
///Recognizement System. Fixed a bug into BASECALC_SOMMAPRIMINNUMERI operation, which didn't display correctly    ///
/// the Sum Result. Enhanced RANDOM-SEED Initializing and Management System, by redefining its default init       ///
/// type into time_t, because actually, starting_random_seed captures the Current Time by using time(NULL) time.h ///
///                          defined function. Enhanced INPUT/OUTPUT Formatting System.                           ///
///                       Redefined Some Program Structures Fields, by using Bit-Fields                           ///
/// assignments. This will not increase in considerable way the Program Speed, but it is a good approach If I     ///
/// i will introduce new array/tables management systems, whose raws and columns Dimensions ID will be probably   ///
///                              incapsulated into BitFields-ful structures.                                      ///
/// Eliminated AUTO_TRANSPOSE Setting, used in the first versions of the Suite for debugging and speed-testing    ///
///  MATRIX_POWER MACRO subProgram. Now it's possible to execute matrixProduct() function and all of its Macro    ///
///                                          with incredible Speed!!                                              ///
///                                      CHANGELOG v5.0 (06/08/2013)                                              ///
/// Totally optimized code and Fixed every imaginable bug into the program, at least the ones that afflicted      ///
/// the program since this update. Added a Fabulous C Parser System and inline functions solver. That's EXPREVAL. ///
/// Now you can manage, calculate various EXPR Expressions and store their results into run-time variables!!!     ///
/// Yeah this isn't a joke. You can perform a Variable List, simply by storing multiple variables with different  ///
/// textual identifiers. And that's not all. You can also store this Varlists into particular files called .msvl, ///
/// (which means Math Suite VARLISTS). It was added an entire program section to manage these MS Environments,that///
/// is an alias for Variable Lists. Variable saving, mantaining and checking processes is totally provided with   ///
/// many hard-coded scripts which, by following the Program Settings Instructions related to this system, store   ///
/// them into the respective files on each expr calculations, or let it onto memory until program is requested to ///
/// exit and to flush the informations on the Disk. Another system of this genre was implemented for Matrix File  ///
/// Extraction Feature. Right as the Varlists, you can store matrixes into Files and later reading them from the  ///
/// Disk. Another program section, located into Linear Algebric Operations, was added to manage Matrixes.         ///
/// Each Matrix printed is saved in memory. You can equal it to the Current Matrix, simply by typing "get" in     ///
/// each ExprEval INPUT Requesting phase. For using the Current Matrix on a calculation which requires a Matrix,  ///
///      simply do type "set" and Program will equal the Working Matrix of the related subprogram to the C.M.     ///
/// However, you have to stay alert because some Programs, you know, requests only square or nXm strong-checked   ///
///     Matrixes, and if you try to equal the W.M. to a C.M. which doesn't have correct requirements for that     ///
/// subprogram, its execution will promptly went discarded and an Error Message will be shown to you. Added also  ///
/// a Log Tracking system, which allows you to register every text-line shown in the program at run-time. Various ///
/// Log files could be loaded in memory, as the previous Lists Item Types upwriten, like Environments and Matrixes///
/// .There is a Current USRLOG File in which the User Executions (like Expr/Matrix calculations will be stored)   ///
/// and also different .LOG Files that are ready to swap with the Current one if you desire it, in any moment.    ///
/// There is also a SYSLOG Management System that stores in a selected .LOG File various System Executions and    ///
/// Settings Changings instead. Default SYSLOG name is (syslog.LOG). However, at the First Execution of the suite,///
/// it is not activated. To enable its execution and functionalities, you must confirm it at the Relative SYSLOG  ///
/// Management Program Section, where also you can set its bufferlen, edit it, rename it and also swap it with    ///
/// with others .LOG Files saved on Disk. There is also a third Lists Item Type, the LAYOUTS one. This system will///
/// allow you to stores Settings List into an .INI File, exactly with the same structure of config.INI one. In the///
/// related Program Section, located into Edit Program Settings, you can Select the CURRENT SETTINGS LAYOUT, or   ///
/// do perform the normal C.R.U.D. (Create, Read, Update and Delete) operations on different Layout previously    ///
/// loaded into the memory. Lists System Items are stored in Heap Memory region, and their storage method is based///
/// on the normal, well-known Double-LINKED-LISTS System, that links each nodelist to the next and prev ones.     ///
/// However, this Program is equipped with a Fast-Recognizement System which searches the Element you want by     ///
///    starting Search in each direction and breaks it when the Items is found. And Program also is based on      ///
/// Fast-Intelligent-ItemID indexing system, which will returns you the Structure contained the pointer to the    ///
/// requested Data, by taking the shortest route in memory. Introduced a new Settings Memory Storage System, which///
/// is based on Bitmasks Improved Usages to abstract Bitflags system one. So a single unsigned variable, obviously///
/// incapsulated into the suite structure is used instead of n-settings bool variables, which will exagerately use///
/// the stack memory. You can also now decide if you want to storage suite main structure into stack or heap,     ///
/// respectively by managing the STACKALLOC Object Macro Definition. Heap Alloc hasn't been tested yet, but you'll///
/// gain a considerable speed increase in theory, at the expense of the Numbers of Lists Items loadable in Heap   ///
/// memory, such as ENVS, MATRIXES, LOGS and LAYOUTS. Added a new system (WINDOWS Only) of Colors Management.     ///
///     You can load an .INI file containing informations about the use of the available Colors in the program.   ///
/// You can view colors.INI default Colors File in order to know how to use this Feature. Added also a new system ///
/// that allows you to load different Lists Items (ENVS, MATRIXES, LOGS and LAYOUTS), simply by putting its fnames///
/// in a .msinf file. Program will load the Items Instances after recognizing them by its extensions. Default ones///
/// I remember, are .msvl for VARLISTS, .msmat for MATRIXES, .LOG for Logs, .INI for LAYOUTS. An example-default  ///
/// MathSuite Information File is autorun.msinf, located in the main root of the program.                         ///
///      Added a very complex Matrix Back-Tracking system, that allows you to back-track to a previous raw/column ///
/// simply by typing "back", while inserting a new element into the Matrix and staying in ExprEval Parsing System ///
/// Settings ENABLED. Code Optimization now includes the exclusive use of pre-increment operator instead of the   ///
/// post-increment one, because this will remove the very-expensive caching mechanism of variable previous value. ///
/// Optimized also all loops: reduced and eventually deleted overhead in the 90% of the programs loops, by using  ///
///         functions pointers instead of performing every time an 'if' control structure et similia...           ///
/// Hard-Prototypized every programs functions, also the private static void ones, by adding attribute lists on   ///
/// each of their declarations. Now it's easily visible whether a function is a system function (so, essential and///
/// very important for having program correctly working), simply by checking if it has __system attribute, or if a///
/// function is a static void (private) or simply not much portable (because is strongly hard-coded for MathSuite ///
/// calculus environment by checking if it has _MS__private attribute. For safe-portable functions you must see   ///
///  for the __export attribute, which (untested proposition) causes defined behavior on some particular IDEs...  ///
/// Added a very powerful and memory-leak PROOF Memoizer System, which allows you to memorize results of the main ///
/// optimized functions of MSCEnv (MathSuite Calculus Environment), such as Fattoriale and Fibonacci, one time    ///
/// they are computed. However, the system will work until a certain index, whose value is stored in the Current  ///
/// SETTINGS Layout. After that, every results of a specific index is calculated normally without mem optimization///
///    . Added ADVANCED_CALCULATOR subProgram, which presents these Children Subprograms: SecondGrade Equations   ///
/// Solver, that by inserting a 3x1 Matrix containing a, b, and c coefficients of ax^2+by+c linear equation std   ///
/// form, allows you to Solve this equation, as the name suggests. There are also simple Complexes Sum and Product///
/// subprograms, which functionality is obviously known. Also there is the possibility of Getting a Formatted Date///
/// in function of a numeric one in the standard EU Format: DD/MM/YYYY. Other subprograms are: Newton Difference  ///
///    Tables, Lagrange Interpolation, Greatest Eigen Value, Function Integration (DEFINED Integral (Riemann) and ///
/// somehow functions must be chosen by a default list), Straight Line Fitting, Parabolic Curve Fitting and also a///
/// very useful and powerful Linear Systems Solver, which allows you to solve multiple equations simply with the  ///
/// entering of the respected System Matrix, which consts of (n+1) Columns, where n is the number of simultaneous ///
/// equations that inhabitates the System. Added a CmdLine RSystem, which emulates the Windows NT MS-DOS Based CMD///
/// Prompt. Obviously, command are Math Orientated and reflects totally the Suite Standard Functions. You can     ///
/// insert multiple parameters separated by a space (seperator system and management will be enhanced in future   ///
/// version of the program). To this feature is dedicated an entire section whose access point is located in the  ///
/// Main Menu (MathSuite Scripts Management). You can enter a CommandLine directly by the program and you can also///
/// load MathSuite scriptfiles, whose extension is tipically .mhss or .msscript (which remembers more than other a///
/// Windows .bat Batch File). This file must contain in each line a MSCenv Command. In order to know usages of the///
/// various commands, there is a dedicated subprogram in this SPS Section (Show Command Usages). You can obtain   ///
/// a singular command usage by inserting as a [STRING] param the name of the CMD (cmdname) or also you can see   ///
/// the entire CMDLINE MacroList, with the relative usages. If you fails to write a command after entered it, MS  ///
/// environment will default show you its usage and re-link you at the Standard CMDLine Program. Improved the     ///
/// inlining of some short-sized functions and enhanced the MINMAX MS Native Embedded powerful system. Replaced   ///
/// the deprecated BubbleSort algorythm to ordinate a vector with the qsort C Native Libraries implementation one ///
/// Removed the Descending Order Option in Matrixes Sort subProgram, which I think useless and stacksize-eater.   ///
///  Removed the max_raws and max_columns management Program located into change_settings.c files, cause these    ///
///    values seems to be more significant defined as Object Macros, since It must not be changed at run-time.    ///
/// Changed the type of some derivate functions like fasum, fnnsum, fsum from uint64_t (unsigned long long) to    ///
/// ityp (Standard MS double Type). Added a basilar Printing System for ENVS, MATRIXES, LOGS, SysLOG and LAYOUTS, ///
/// available both for WINDOWS and UNIX systems. The unique difference between them is that whilst in Windows you ///
/// can select the Device, in UNIX you cannot and spooling settings are default-set. However, this feature hasn't ///
/// been tested yet). Enhanced SecurityCheck, now the Stack is high-protected by some controls that are strongly  ///
/// hard-coded into program. The Exit Button is now disabled during Matrix Inserting Processes. This preserves the///
/// unexpected exits from the program when malloc, calloc and free NATIVE malloc.h library functions are working  ///
/// . This system has the advantage also of prevent the memory getting too segmentated (works only on WINDOWS).   ///
/// IMPORTANT ADVICE: Program isn't completely portable now. For getting the Suite right-working on UNIX based    ///
/// System you must apport some changes to the source code and recompile-it under UNIX Environments, by following ///
/// also the Build messages advices that compiler or IDE provides for you... I recommend Code::Blocks for WINDOWS ///
/// users. Now it's hardly recommended to exit by the related Option that Programs offers to you. This will allows///
/// you to safeExiting the program and NULL-Setting all the pointers (could apparently be useless this feature but///
/// it isn't really) and some particular variables that deals with Memory Management (such as the Dynamic Array of///
///                             of Pointers for the MS Embedded Memoizer System).                                 ///
/// Added a new system of calculating the Determinant of a Square Matrix, by upper-triangularizing it with the    ///
/// respective matrixUTConv, which I implementated under checking Bibek Subedi programming-techniques.com various ///
/// pieces of code. This technique has been also used in other Suite SubPrograms, like Linear Systems Solvers     ///
/// (also partially taken from that website). Now Rank Calculator Modules uses Jacobi SVD (Singular Value         ///
///       Decomposition, which I taken from the web and optimized & converted to get working on my program.       ///
/// Substantially, RkCalculator now takes the Filed Singular Values Vector from this function dsvs(...) and count ///
/// its non-zero elements. Now, when using DMA (Dynamic Memory Allocation), it's used calloc function instead of  ///
/// malloc one. This was done for Security Reason, because, even If very-unlikely there could be a bug on Matrix  ///
/// Calculations, It can't damage memory or brings the program to freeze or crash, because a nan value, even If   ///
/// wrong, could block its execution. However, I planned to insert a very powerful saturate addiction system that ///
/// blocks execution before computation. Even if there is a Basilar Overflow Checking System, implemented in      ///
/// programs former versions, and even If it has been also enhanced in this update, could be wrong and so, could  ///
/// inform you too numbers later of the Inconsistent Result you get. So, even If partially deprecated, this base  ///
/// system is efficient and powerful meantime, because a Saturate Addiction (or Algebric Sum somehow) requests    ///
/// two or three control structures that could (Even If in a not considerable way) affect the Execution Time.     ///
/// Added the possibility of viewing the DATE and TIME on every single ProgLine when the related Bool Setting Var ///
/// is Set in the Current Layout. You can also decide to see the Average Time of each ExprEval Computation or the ///
/// Execution Time of Every SubProgram, including the Access Point ones (the ones that allows you to get in a     ///
/// program menu section), simply by enabling respectively the DiffTime and ExecTime Bool Vars ones in Current    ///
/// Layout. Enhanced the Mechanism of Selecting the Item from a constant strings list. This is also used for      ///
/// choising the itemID of a Lists Items when handling __lmp_prog (Lists Manager Program) and means more powerful ///
/// optimizations given to the program. Completely redesigned the Prime_N_Number functions, improving a new       ///
/// mechanism of checking whether a Numer is prime or not. They are available two algorythm to perform the Sprgrm ///
/// , respectively isPrimeHISP and isPrimeLIFP (the first works on High Iterations Low Process and the second     ///
/// on Low Iterations, Fast Process). When Lazy-EXECUTION CurrentLayout Bool Var is Set, isPrimeLIFP is used,     ///
/// otherwise HISP method will be. Enhanced the Presentation Layer and worked so hard to get a primitive and very ///
/// confidential INPUT/OUTPUT Interaction. In order to increase the speed when inserting a new Matrix, ExprEval   ///
/// Parser System combined to my work gave you the possibility to use the different Environments Variables also   ///
/// when performing Automatic Matrix Inserting!!! This means that you can calculate for example, three variables  ///
/// with BASECALC subProgram and use their values with other computations directly performed into the [2,3] item  ///
/// of a Matrix getting run-time inserted. Even if is strongly recommended to use ExprEval to do computations     ///
/// involving complex EXPR, the old BASECALC method is still present. If you want to use It, simply disable       ///
/// the Parsing for BaseCalc Program Setting Option of the C.L. (Current Layout). Basecalc Traditional way of     ///
/// calculating isn't really useless, because the OVERFLOW_CHECKER and DOMAIN_CHECKING systems, works only in the ///
///   traditional mode. The Math Errors routine with the Parser ON are managed directly by ExprEval, so by its    ///
/// errors manager own .c files. Somehow you must not get yourself involved with Math Errors very frequently,     ///
/// when a 'nan' value is displayed (some bug for example), means that you overflowed an ityp value (still EX..)  ///
/// Introduced also a simple ALIASES System that allows you to have more Aliases for the same Traditional BCALC   ///
/// command. This system, for simplicity reasons, is suitable only with having the source code re-compiled. All   ///
/// these features that require re-compiling the program are however not useless since MathSuite Environment is   ///
/// aimed to Open-Source Programming, Portability and WIN-UNIX oriented. Added a basilar but powerful catchPause  ///
/// system, consisting of giving user the possibility of blocking a subProgram execution while for example,       ///
/// printing a vector or File content, simply by pressing CTRL+C. This could be useful when listing prime numbers ///
/// whose execution and numbers printing could easily overflows your stdout buffer too early and not more giving  ///
/// you possibility to see previous results. Now you can Flush also the stderr file buffer. Added three new sprogs///
/// into algebra.c subfile: SMatrix LU Decomposition, Matrix Ill Condition Checkin, Square Matrix Norm Calculation///
///                                 and also SVD Matrix Decomposition Vector Viewer                               ///
/// 				Syntax optimized and Code optimized [06/08/2013].				  ///
///                                         CHANGELOG v5.2 (05/08/2014)					 	  ///
/// Fixed a lot of little bugs that afflicted the program in some non-standard situations. Syntax optimized and   ///
/// Code highly optimized. checkItemTypeByExtension has been redesigned in order to increase speed execution      ///
/// while checking items extensions, for example in loading a .msinf startup file. Fixed a bug in the Primality TS///
/// Engine, by adding the "2" Prime Number, considered obvious in earlier versions. Now the default mode of suite ///
/// main vars structures' allocations is Heap Mode. (STACKALLOC Macro isn't defined by default). Now the program  ///
/// seems to be more efficient and speedy, but probably the max logical number of items to allocate in memory has ///
/// been reduced, due to the main structure allocation. However, with the increasing of the executable dimension  ///
/// this has been inevitable. Added new 105 functions, respectively 95 trigonometric and pseudo-trigonometric     ///
/// (hyperbolic), and 10 exponential and logarithmic. These has been added both in Base Calculator and in MSCenv  ///
/// Parsing System. These are: hsin(a), hsinh(a), qsin(a), qsinh(a), hcos(a), hcosh(a), qcos(a), qcosh(a), hsec(a)///
///   hsech(a), qsec(a), qsech(a), hcsc(a), hcsch(a), qcsc(a), qcsch(a), htan(a), htanh(a), qtan(a), qtanh(a),    ///
///   hcot(a), hcoth(a), qcot(a), qcoth(a), vsin(a), vsinh(a), cvsin(a), cvsinh(a), hvsin(a), hvsinh(a), qvsin(a),///
///   qvsinh(a), hcvsin(a), hcvsinh(a), qcvsin(a), qcvsinh(a), vcos(a), vcosh(a), cvcos(a), cvcosh(a), hvcos(a),  ///
///     hvcosh(a), qvcos(a), qvcosh(a), hcvcos(a), hcvcosh(a), qcvcos(a), qcvcosh(a), esec(a), esech(a), ecsc(a), ///
/// ecsch(a), hesec(a), hesech(a), hecsc(a), hecsch(a), qsec(a), qsech(a), qcsc(a), qcsch(a), sinc(a), sinch(a),  ///
/// hsinc(a), hsinch(a), qsinc(a), qsinch(a), cosc(a), cosch(a), hcosc(a), hcosch(a), qcosc(a), qcosch(a), secc(a)///
/// secch(a), hsecc(a), hsecch(a), qsecc(a), qsecch(a), cscc(a), cscch(a), hcscc(a), hcscch(a), qcscc(a),qcscch(a)///
/// tanc(a), tanch(a), htanc(a), htanch(a), qtanc(a), qtanch(a), cotc(a), cotch(a), hcotc(a), hcotch(a), qcotc(a),///
///     qcotch(a), logc(a), log10c(a), log2c(a), exp10(a), expc(a), exp10c(a), exp2c(a), log1p(a), log1pc(a).     ///
/// These functions have been introduced also in the PRELoaded Functions Integration System, so giving the        ///
/// 	to integrate them, by using the available integration methods: the Simpson and the Trapezoidal.	  	  ///
///                    Corrected some translation error in environment printf phrases.                            ///
///                                     CHANGELOG v5.21 (06/08/2014)                                              ///
/// Added the Kronecker Product feature in the Algebra Operations Group. Fixed a bug in the main code, that       ///
/// have been crashed the program if it received a one-length string as argv[1]. Fixed a bug in the Matrix Sum    ///
/// subProgram, that maybe could cause an instant crash if the Second Matrix Entering process failed. (U.T.)      ///
/// Renamed some descriptions file in the apposite folder, that caused a bug if the program requested descriptions///
///                         for the Matrices Product and Matrices Sum subPrograms.                                ///
///                  Fixed some bugs and code optimized (these are minor things, however...)                      ///
///                                     CHANGELOG v5.50 (16/08/2014)                                              ///
/// Fixed and Optimized code. Fixed a bug into the Kronecker Product, that in some random cases could crash the   ///
/// program because an heap overflow at the moment of equaling the Last Matrix Printed (I remember that Kronecker ///
/// Product is slightly heavy for Heap Memory and also for CPU (it takes 4 nested for cycles to work)). Introduced///
/// Matrix Kronecker Power feature, that allows you to Kronecker-raise the inserted matrix to the inserted power  ///
/// exponent, like in the normal Matrix Power subProgram. Introduced the semi-factorial function. Its name is:    ///
/// "sfact", so I renamed correctly the Stabilizer Factor Function into "stabfact" or similar (see its ALIAS      ///
/// MACROS for further informations. This function works upon Factorial Memoizer System, but it hasn't an own     ///
/// Memoizer Engine, because its non proper standard complexity. Maybe i will introduce later a new MemoizerEngine///
///     for this purpose. Added the related function "sfasum". Fixed a bug into the appendTimeToString function,  ///
///   now instead of overflowing the stdout buffer and making illegible the Current Textual Environment, it will  ///
///append Current Time to each textline only if the Line Color if different from the previous one. But in some    ///
/// rare cases, like in MathSuite's logo printing process it isn't so useful. Automatized some Heap Dynamic Memory///
///       Allocation Management Processes, by incorporating some default routines into Matrices IO Functions      ///
///(for example matrixAlloc into insertNMMatrix). Fixed a bug (indeed it isn't a proper bug but my old            ///
/// technical choice) that have not included the number "2" into Prime Numbers Selection or Primality Tests.      ///
/// Introduced the Simplex Method, precisely the Non-Dual Simplex Method for solving PL Problems with non-negative///
/// variables. You can solve a simple PL problem with these characteristics, simply by inserting a Matrix         ///
/// containing Constraints Coefficients, and let the last raw contains the Functions Coefficients indeed, followed///
/// by an unique 0 element to align Matrix Dimension. Successively you have to inform the program about the       ///
/// Constraints Types, simply by inserting a n-dimensional Vector (following Program Instructions is relatively   ///
/// simple), whose i'th element has to be 0 if the i'th constraint is in the type of '<=' condition. It has to be ///
///                                     an integer different from 0 otherwise.                                    ///
///Added the new Algebra Selection system, that allows you to operate in different                                ///
/// algebras while performing Linear Algebra Operations such as Matrix Sum, Matrix Product or Matrix Power, Matrix///
/// Kronecker Product or Matrix Kronecker Power. You can choice from Real Numbers (naturally, the default option) ///
///, Complex Numbers, Quaternions, Octonions and Sedenions!!! Yeah you've heard me really!! You can also perform, ///
/// for example, a Kronecker Product with two Matrices whose elements are Sedenions! And naturally, the program   ///
/// will perform the atomic algebric operations between elements by following Selected Algebra Rules. But to do   ///
/// this, you have to insert: one Matrix in any case containing all the Real Parts of the numbers, and many       ///
/// matrices as many the Imaginary Parts of the Current Algebra are, simply by filling the respective matrices    ///
/// with the Imaginary Parts Coefficients of the numbers. The Imaginary Parts are: 1 for Complex Numbers (i), 3   ///
/// for Quaternions (i, j, k), 7 for Octonions (e1, e2, e3, e4, e5, e6, e7) and 15 for Sedenions (e1, e2, e3, e4, ///
/// e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15). Note that you can also perform single atomic operations    ///
/// between complex and hyper-complex numbers in the homonymous subProgram, that I've updated in this release, in ///
///           the Advanced Calculator Section, the same one in which Simplex Method new feature is located.       ///
/// Introduced the Tensor Sum subProgram, that allows you to sum two n-dimensional Tensors by performing n sum    ///
/// between the respective Component Matrices. Fixed a bug into Current Layout Parser, that because an integer    ///
/// overflow, program didn't read some Layout Values correctly, such as Matrix Max Raws, Matrix Max Columns, etc..///
///                                     CHANGELOG v5.55 (20/08/2014)                                              ///
/// Code optimized and Fixed some minor bugs. Fixed a bug into printMatrix function, that didn't allow you to     ///
/// equal the Snapshot Matrix (lmpMatrix) to the effective Last Matrix Printed, and even in some cases it might   ///
/// have crashed the program (for example in the Matrix Kronecker ORIGINAL Release (v5.21)). Fixed a bug into the ///
/// matrixTransposition function, that crashed the program If someone attempts to insert correctly the Matrix.    ///
/// Introduced Multi-Threading and Extensive Multi-Threading by using OpenMP (Open-MultiProcessing) libraries and ///
/// functions. It has been possible to parallelize all the critical sections which contain several independent    ///
/// iterations. Added the possibility to switch between Multi-Threading and Extensive Multi-Threading: the first  ///
/// one system, forces in some cases the creation of N threads; but if this N number exceeds the number of an     ///
/// 		hard-coded constant, defined both in defaults.h and successively in dutils.h, named		  ///
/// MIN_EXTENSIVE_MULTITHREADING_CORESNO, then let OpenMP manage the threads creation process. Elsewhere, with the///
/// 	Extensive Multi-Threading system enabled, program could appreciate the parallel power of more than, for   ///
/// example, 4 cores, because in this case, OpenMP will be forced to create exactly N threads. This operation is  ///
/// very important, for example, when operating in Different Algebras whilst performing Linear Algebra operations,///
/// because Algebra Units are all powers of 2 (example 2,4,8,16), exactly like the number of cores of a modern PC.///
/// The SemiFactorial's function implemention has switched to an iterated algorythm. Moreover, they have been     ///
/// created two different buffering channels, that serve as the new Memoizer Engine for SemiFactorial function.   ///
/// The first pipe deals with Even number, the second pipe with Odd numbers, because Memoizer System is strongly  ///
/// index-based and the SemiFactorial function itself needs two different type of previous values for computing   ///
/// 	its next value: respectively to calculate the semifactorial of a 2n typed number, it needs recursively	  ///
/// numbers in the form of (2n-2)!!, which is: 2(n-1)!! that are surely Even Numbers. Elsewhere, if n is an Odd   ///
/// 			number, it will require only Odd previous numbers to calculate its return value.	  ///
/// Introduced Stirling's Approximation in the Factorial function implementation: if the number of which we want  ///
/// 	 to calculate the Factorial, exceeds a runtime constant located into the Current Layout and named:	  ///
/// min_stirlingrequires_number, then fact(n) function will call the stirling(n) one, that is supposed to reduce  ///
/// considerably the Overflow Error. This optimization is valid also for the SemiFactorial function, but only for ///
///    Even number, due to the existant relation between Factorial and SemiFactorial: (2^n)*n! = (2n)!! ==>	  ///
/// 	    	(2n)!! ---> (2^n)*stirling(n). This is possible if 2n > 2*min_stirlingrequires_number	          ///
/// Renewed all the Selection Processes into Settings Manager: I replaced the deprecated do ... while method  	  ///
/// containing toupper(getch()) in its body, with the more powerful and secure selectListItem(...) function. 	  ///
/// This change has been applied to: "Change Algebra", "Empty Memoizers Buffers" and "Empty Buffers" subPrograms. ///
///                                         Code totally Optimized.                                               ///
///                                     CHANGELOG v5.60 (23/08/2014)                                              ///
///   Fixed some minor bugs and Code totally optimized. Introduced some statistic functions, such as: Variance,   ///
/// Covariance, Standard Deviation, Outlier Test, First Quartile, Third Quartile, respectively with these 	  ///
/// identifiers: var, cov, stddev, otlr, otlr2, q1, q3. Fixed a bug in the Mediana Function, which wasn't 	  ///
///doing correctly its job, because an inverted test of a compact boolean expression. Parallelized some code areas///
///	Optimized some iterations structures in the Base Calculator, precisely in the Statistic V.A.R.		  ///
/// (Vector-Accumulating-Requires) Functions, like media, etc... Removed all the Temperature Conversions Functions///
/// and also the Speed Conversion Function, because even If they're pretty smart and useful function, they're also///
///					too specific for a Math Related program.				  ///
///                                     CHANGELOG v5.70 (29/08/2014)                                              ///
/// Discarded the .INI system for gathering informations about Program Settings and Colors Settings, and replaced ///
/// with the XML metalanguage for the same purpose. The Layouts System works exactly as the previous versions.    ///
/// Renamed all the extensions of the program: .msvl to .vlf (it stands for Variables List File), .msmat to .mf   ///
///   (Matrices Files), .mhss to .mss (MathSuite Scripts), .msinf to .lf (List Files).				  ///
/// Now the presence of autorun.lf startup File isn't requested anymore, because even If the program is not be    ///
///    able to load at least a Layout File, then It will create a default settings.xml or colors.xml file. 	  ///
/// Now it is possible to save the Exit Char and the Outlier Constant into a Layout File. Added the following	  ///
/// commands: ec, oc, map; respectively the first one changes Exit Char, the second will do the same for the      ///
/// Outlier Constant and the third command, if Base Calculator Parser is enabled then It will apply a function to ///
/// an inserted value. The function is identified by a constant inserted by the user, named FID; if Basecalc	  ///
/// Parser is enabled otherwise, then it will apply to the function to a vector. This vector is inserted by       ///
/// following the V.A.R. functions standard method (by accumulating a vector into stack and then free it when you ///
/// have inserted all the elements). When you're done, then you have to insert a different value of the FID in    ///
/// order to stop Vector Accumulating Process. Otherwise you have to insert the same initial value of FID to      ///
/// continue. Fixed some minor bugs and Code highly optimized and stabilized the program run-time execution.	  ///
///                                     CHANGELOG v5.75 (01/09/2014)                                              ///
/// Added Permutations, Permutations with Repetitions, K-Permutations, K-Permutations with Repetitions and        ///
/// Combinations with Repetitions, with the respective commands: "perm", "permrep", "kperm", "kpermrep", "combrep"///
/// 			Removed "Vector per Vector" MACRO SubProgram, due to an undefined bug.			  ///
///                                     CHANGELOG v6.00 (04/09/2014)                                              ///
/// Now the Matrix Base Type with double-reference pointer of double base-type has been changed to a single-ref   ///
/// pointer. Now the maximum reference depth level present in my program is three (for the Tensors Sum feature).  ///
/// Optimized program and deleted all metadata management system related to the previous method of allocating a   ///
///  new chunk of memory every time that a matrix needed a new raw. Added the "restrict" qualifier type to the    ///
/// single-reference M.B.T. (Matrix Base Type) in the first-level functions arguments. Fixed a bug into the       ///
/// Complex and HyperComplex both Sum and Product subPrograms, that didn't allow to view correctly the Result of  ///
/// the Operation due to an array field accessing problem. Now the S.E. (STABILIZER+INITIALIZER) System, known    ///
/// also as M.S.D.M.A (mathSuite Dynamic Matrix Allocation) has been deprecated due to some problems occurred with///
/// 	the introduction of the Single Reference M.B.T.. So, the Stabilizer System works only with vectors.       ///
/// Added "qsum" and "qprod" inline functions, that performs respectively Sum and Product between two Quaternions.///
/// 				    Code Totally optimized and fixed also some minor bugs. 			  ///
///                                     CHANGELOG v6.50 (10/09/2014)                                              ///
/// Fixed a bug into the extensive multi-threaded matrix Sedenions Product routine, which didn't perform correctly///
/// the operation on the e6 base-field of the matrix. Introduced new 154 complex functions, respectively 120      ///
/// trigonometric and pseudo-trigonometric (hyperbolic), 23 exponential and logarithmic, 5 about complex          ///
///   argument and complex and hypercomplex numbers absolute values, 2 are quaternions inline Addition and        ///
/// Multiplication and the last 4 are the extensions of the 1P Logarithmic Functions (inclusive of the Cardinal   ///
/// Versions) at the other bases. All these functions are respectively: chsin(a,b,&rrp,&rip),chsinh(a,b,&rrp,&rip)///
/// cqsin(a,b,&rrp,&rip), cqsinh(a,b,&rrp,&rip), chcos(a,b,&rrp,&rip), chcosh(a,b,&rrp,&rip), cqcos(a,b,&rrp,&rip)///
///cqcosh(a,b,&rrp,&rip), chsec(a,b,&rrp,&rip), chsech(a,b,&rrp,&rip), cqsec(a,b,&rrp,&rip), cqsech(a,b,&rrp,&rip)///
/// chcsc(a,b,&rrp,&rip), chcsch(a,b,&rrp,&rip),cqcsc(a,b,&rrp,&rip), cqcsch(a,b,&rrp,&rip), chtan(a,b,&rrp,&rip),///
///chtanh(a,b,&rrp,&rip), cqtan(a,b,&rrp,&rip), cqtanh(a,b,&rrp,&rip), chcot(a,b,&rrp,&rip), chcoth(a,b,&rrp,&rip)///
///       ,cqcot(a,b,&rrp,&rip), cqcoth(a,b,&rrp,&rip), cpxvsin(a,b,&rrp,&rip), cpxvsinh(a,b,&rrp,&rip),          ///
///	   ccvsin(a,b,&rrp,&rip), ccvsinh(a,b,&rrp,&rip), chvsin(a,b,&rrp,&rip), chvsinh(a,b,&rrp,&rip),          ///
///      cqvsin(a,b,&rrp,&rip), cqvsinh(a,b,&rrp,&rip), chcvsin(a,b,&rrp,&rip), chcvsinh(a,b,&rrp,&rip),          ///
/// 	cqcvsin(a,b,&rrp,&rip), cqcvsinh(a,b,&rrp,&rip), cpxvcos(a,b,&rrp,&rip), cpxvcosh(a,b,&rrp,&rip),         ///
/// 	ccvcos(a,b,&rrp,&rip), ccvcosh(a,b,&rrp,&rip), chvcos(a,b,&rrp,&rip), chvcosh(a,b,&rrp,&rip),	    	  ///
///	cqvcos(a,b,&rrp,&rip), cqvcosh(a,b,&rrp,&rip),chcvcos(a,b,&rrp,&rip),chcvcosh(a,b,&rrp,&rip),	    	  ///
///cqcvcos(a,b,&rrp,&rip),cqcvcosh(a,b,&rrp,&rip),cesec(a,b,&rrp,&rip),cesech(a,b,&rrp,&rip),cecsc(a,b,&rrp,&rip),///
///	  cecsch(a,b,&rrp,&rip), chesec(a,b,&rrp,&rip), chesech(a,b,&rrp,&rip), checsc(a,b,&rrp,&rip), 		  ///
///checsch(a,b,&rrp,&rip), cqsec(a,b,&rrp,&rip), cqsech(a,b,&rrp,&rip), cqcsc(a,b,&rrp,&rip),cqcsch(a,b,&rrp,&rip)///
///csinc(a,b,&rrp,&rip), csinch(a,b,&rrp,&rip), chsinc(a,b,&rrp,&rip),chsinch(a,b,&rrp,&rip),cqsinc(a,b,&rrp,&rip)///
///	  cqsinch(a,b,&rrp,&rip), ccosc(a,b,&rrp,&rip), ccosch(a,b,&rrp,&rip), chcosc(a,b,&rrp,&rip),		  ///
///	  chcosch(a,b,&rrp,&rip), cqcosc(a,b,&rrp,&rip), cqcosch(a,b,&rrp,&rip), csecc(a,b,&rrp,&rip),		  ///
///		csecch(a,b,&rrp,&rip), chsecc(a,b,&rrp,&rip), chsecch(a,b,&rrp,&rip), cqsecc(a,b,&rrp,&rip),      ///
///	     cqsecch(a,b,&rrp,&rip), ccscc(a,b,&rrp,&rip), ccscch(a,b,&rrp,&rip), chcscc(a,b,&rrp,&rip),    	  ///
///	   chcscch(a,b,&rrp,&rip), cqcscc(a,b,&rrp,&rip),cqcscch(a,b,&rrp,&rip), ctanc(a,b,&rrp,&rip),       	  /// 
///	    ctanch(a,b,&rrp,&rip), chtanc(a,b,&rrp,&rip), chtanch(a,b,&rrp,&rip), cqtanc(a,b,&rrp,&rip),	  ///
///cqtanch(a,b,&rrp,&rip), ccotc(a,b,&rrp,&rip),ccotch(a,b,&rrp,&rip),chcotc(a,b,&rrp,&rip),chcotch(a,b,&rrp,&rip)///
///cqcotc(a,b,&rrp,&rip), cqcotch(a,b,&rrp,&rip),clog(a,b,&rrp,&rip), clogc(a,b,&rrp,&rip), clog10(a,b,&rrp,&rip),///
///clog10c(a,b,&rrp,&rip), clog2(a,b,&rrp,&rip),clog2c(a,b,&rrp,&rip), cexp(a,b,&rrp,&rip), cexpc(a,b,&rrp,&rip), ///
///cexp10(a,b,&rrp,&rip), cexp10c(a,b,&rrp,&rip), cexp2(a,b,&rrp,&rip),cexp2c(a,b,&rrp,&rip),clog1p(a,b,&rrp,&rip)///
///	clog1pc(a,b,&rrp,&rip), clog101p(a,b,&rrp,&rip), clog101pc(a,b,&rrp,&rip), clog21p(a,b,&rrp,&rip),        ///
///   clog21pc(a,b,&rrp,&rip), carg(a,b,&rrp,&rip), cabs(a,b,&rrp,&rip), qabs(a,b,c,d), oabs(a,b,c,d,e,f,g,h),    ///
///       sabs(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p), qadd(a1,b1,c1,d1,a2,b2,c2,d2,&rrp,&rip,&rjp,&rkp),	 	  ///
/// 	qmul(a1,b1,c1,d1,a2,b2,c2,d2,&rrp,&rip,&rjp,&rkp), log101p(a), log101pc(a), log21p(a), log21pc(a).	  ///
///    The upwriten functions perform calculations and store the Results into the ampersanded variables.	  ///
/// Fixed a bug into "cadd" and "cmul" commands, that didn't allow to store the correct result into the rightmost ///
/// variables arguments. Fixed a bug into the Complex and HyperComplex respectively both Addition and 		  ///
/// Multiplication subPrograms, that didn't allow to format correctly the results. Totally rewritten the Domain   ///
/// Checking System related functions and so enhanced precision in all the trigonometric functions. Now even the  ///
/// Basic Calculator have the new Complex Trigonometric, Exponential, Logarithmic and Miscellaneous functions     ///
/// introduced in this version, and the Domain Checking System (D.C.S.) works perfectly also in this Environment! ///
/// Furthermore, it wasn't possible to accorp one complex trigonometric function with its hyperbolic version one  ///
/// 	   one. So, bcilfs.h Basic Calculator HEADER got a bit crowded, but it isn't really a problem. 		  ///
/// Introduced a simple performant OSMM Algorithm for the Square Matrices Multiplication (O.S.M.M. means Optimized///
/// Square Matrices Multiplication) with dimension greater than a constant defined into the settings.xml file     ///
/// (or any valid XML Settings File). Hoever it is possible to change this constant by using its subprogram into  ///
/// "Edit Program SETTINGS" section. This algorithm works also with a BLOCK_SIZE constant that depends on some    ///
/// hardware cache-related specifications. It is also possible to change its value both in XML Settings File and  ///
/// by its inline setting command. Introduced also the Strassen Algorithm for Square Matrices Multiplication      ///
/// with dimension of the form of: 2^n, with n integer. If the "Strassen Optimization" Bool Settings is enabled,  ///
///then the program will use this algorithm instead of the naive classic one whenever the dimension of the actual ///
/// Matrices is greater than a constant defined into the XML Settings File, just as the OSMM Algorithm one. Note  ///
///that if these conditions are true and the dimension of the matrices is greather than the Min OSMM Dimension    ///
/// 	     constant, Strassen has the priority upon OSMM algorithm for square matrices. 			  ///
///					Fixed minor bugs and Code highly optimized.				  ///
///                                     CHANGELOG v6.60 (15/09/2014)                                              ///
/// Introduced new checks into XML functions. Now if during a writing/parsing operation one field of both two XML ///
/// files (settings.xml or colors.xml), is missing, then the program will continue its execution and will show    ///
/// an error message, or eventually, will also set a default value to the interested setting. Fixed a bug into    ///
/// files opening process interface, precisely in Item Selection by Path Modality, that didn't show an error      ///
/// message if the inserted path referred to a non-existent file, or it didn't exit even if you entered the       ///
/// special Exit Char for this mansions ("."). Some outputs has been tightly decorated, though not in an excessive///
/// way. NOTICE: Since this release, the project is completely PORTABLE! It has been compiled under Linux system. ///
///                                     CHANGELOG v6.70 (16/09/2014)                                              ///
/// Replaced the deprecated subprogram "Greatest Eigen Value" of Advanced Calculator section, with the new more   ///
/// performant "Matrix Eigen Values", now placed into Linear Algebra Operations section. This program, as the name///
/// suggests, calculates Eigen Values and Eigen Vectors of a given square matrix. Fixed some minor bugs and C.O.  ///
///                                     CHANGELOG v6.80 (25/09/2014)                                              ///
/// Fixed a bug into the productory(...) ext_math.c function, whose ExprEval inline identifier is "product", which///
/// didn't calculate correctly the result due to a bad initialization of a counter variable. Now this function and///
/// all the functions which was using it, work perfectly. Eliminated the old User Interface related to the Basic  ///
/// Calculator with the Parser Modality disabled. So, reduced the executable size of the 41.5%. Introduced a new  ///
/// 	   statistical function, the "mode", and the code optimized in some critical parts.			  ///
///                                     CHANGELOG v6.85 (04/10/2014)                                              ///
/// Re-written some critical parts of the Lists Items Manager Engine, by using static arrays of functions or vals ///
/// instead of using long cascades of switches. Introduced the CoFactor Matrix subProgram and the related Adjoint ///
/// program, which, as the names suggest, calculate respectively the CoFactor and the Adjount Matrices of a given ///
/// 						  real square matrix.						  ///
///                                     CHANGELOG v7.00 (15/11/2014)                                              ///
/// Introduced the possibility to calculate both the Routh Table and the Jury Table of a n>2 dimensioned RowMatrix///
/// 	 whose elements are the coefficients of the polynom to which you want to apply those Criterions.	  ///
/// 					Fixed some bugs and Code optimized.					  ///
///    This is not a Final Built or Official Release, but (UPDATE) probably this could be the last release I write///
///              in C. A future release could be directly mathSMART Mobile Math Environment.                      ///
///                     You can contact me at: marco_chiarelli@yahoo.it or on the secundary mail:                 ///
/// marcochiarelli.nextgenlab@gmail.com in order to report a bug or simply for sending me an advice that could be ///
///                        useful or could improve the speed or optimize my MSCenv System.                        ///
///  --- This is not a Final Built or Official Release, even If the program seems to be efficient and completely  ///
///                  beta-tested. You can contact me at: marco_chiarelli@yahoo.it or on the secundary mail:       ///
/// marcochiarelli.nextgenlab@gmail.com in order to report a bug or simply for sending me an advice that could be ///
///                        useful or could improve the speed or optimize my MSCenv System.                        ///
//!-------------------------------------------------------------------------------------------------------------- ///
///                                    TO-DO in EVENTUAL NEXT VERSIONS                                            ///
/// Predisposed the project to Join a new powerful Calculus System, that processes an operations by dividing it   ///
/// in specified nibbles and processing them singularly by different processes, that will communicate each them   ///
/// by PID signals sending/receiving. The nibbles quantity depends from the number size/overflow error happening  ///
/// probability.                                                                                                  ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!______________________________________________________________________________________________________________ ///
//!-------------------------------------------------------------------------------------------------------------- ///
///                contact me at: marco_chiarelli@yahoo.it or marcochiarelli.nextgenlab@gmail.com                 ///
///        I'll be glad to fix your scripts or simply to take away your doubts about the program                  ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!-------------------------------------------------------------------------------------------------------------- ///
/// Thanks to giggikr: http://forum.html.it/forum/showthread/t-1374455.html for his own function, cambiabase,     ///
///                         which I renamed, modified and adapted to this program. Thanks to:                     ///
/// http://elite.polito.it/files/courses/12BHD/progr/Esercizi-C-v2_01.pdf for some of their scripts.              ///
/// Thanks to Bibek Subedi, for his invertMatrix function, which I renamed, modified and adapted to this program. ///
/// Link Source: http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html    ///
/// Thanks to Paul Bourke:https://www.cs.rochester.edu/u/brown/Crypto/assts/projects/adj.html for his CoFactor fnc///
/// 	   		     which I modified and adapted to this program.					  ///
/// Thanks to W. Cochran  wcochran@vancouver.wsu.edu for his Strassen Algorithm Implementation, which I renamed,  ///
/// adapted and modified to this program. Thanks also to: Computer Science Division | EECS at UC Berkeley for     ///
/// some notions about Matrix Multiplication Optimizations Techniques: www.cs.berkeley.edu/~knight/cs267/hw1.html ///
/// 	Thanks to: http://cap-lore.com/MathPhys/eigen/j.c for the actual _matrixEigenValues function.		  ///
/// Massive thanks to Brian Allen Vanderburg II for his fabulous C parser and inline functions solver, EXPREVAL,  ///
/// which elegantly gave in theory infinite functionalities and potential to my program. That's the project link  ///
/// with Online Documentation: http://expreval.sourceforge.net/ Thanks to: http://www.cprogramming.com/tips/ and  ///
///             http://stackoverflow.com/questions/599365/what-is-your-favorite-c-programming-trick               ///
///                       http://stackoverflow.com/questions/132241/hidden-features-of-c                          ///
///    that are some websites in which I found a lot of useful C tips and tricks, and they were an important      ///
///  checkpoint for resources retrieving in order to speed-up and optimize my program. Still greatly thanks to    ///
///    Bibek Subedi for his website: http://www.programming-techniques.com/ which put in front of my eyes a new   ///
/// world of C programming. I also recently renewed the program code by improving a lot of his C tricks and tips. ///
/// For example, the upper-triangular Matrixes conversion, which was useful to enhance some functions like det(), ///
///      				sgeqsolver ExprEval inline command, etc.	  			  ///
///	  		Greatly thanks to Daniel Veillard for his fabulous XML Parser, LIBXML2.			  ///
///                 Greatly thanks to Daniel Veillard for his fabulous XML Parser, LIBXML2.                       ///
///  Greatly thanks to vict85 of matematicamente.it Network, for having informed me about the benefits of using   ///
///   generally a single reference for the Matrix Type, like LAPACK and the other Numeric Calculus Environments.  ///
/// 		Thanks to Francesco Palma for reporting me some bugs, and finally, massive thanks to my		  ///
/// Informatic Fundaments Teacher, Mario Alessandro Bochicchio, which gave me a lot of C advices and some general ///
///            tricks and tips, that enlarged my professional informatic horizonts. That's all...                 ///
/*!________________________________________________________________________________________________________________*/
//!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
