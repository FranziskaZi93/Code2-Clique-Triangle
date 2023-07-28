
import robustTwoClub.algorithms.RobustTwoClubAlgorithm;
import robustTwoClub.algorithms.RobustTwoClubAlgorithm.Model;
import robustTwoClub.correctness.CorrectnessChecker;
import robustTwoClub.graph.RtcGraph;


public class RobustTwoClub {

	public static void main (String[] args) {
		
		long start = System.currentTimeMillis();
		
		if (args.length == 0) {
			System.out.println("Called program without arguments. Must at least specify input graph file.");
			System.out.println("For example:   java RobustTwoClub test.graph");
			System.out.println("Input graphs must be either in DIMACS, Metis or clq format.");
			System.out.println("Call 'java RobustTwoClub -help' to list further options.");
			System.exit(0);
		}
		
		if (args[0].equals("-help")) {
			System.out.println("\nThe first argument must always be the name of the input file.");
			System.out.println("\nFollowing options can be added after the filename:");
			System.out.println("-vb        \t Solve for the robust 2-club model of Veremyev/Boginski.");
			System.out.println("           \t (Adjacent vertices also require t neighbors.)");
			System.out.println("-tc        \t Solve for the t-connected 2-club model.");
			System.out.println("           \t (Resulting 2-club has to be t-connected)");
			System.out.println("           \t Note: ");
			System.out.println("           \t If neither the -vb nor -tc flag is set, then the program finds.");
			System.out.println("           \t the largest (t-1)-hereditary 2-club.");
			System.out.println("           \t (A t-robust 2-club is a (t-1)-hereditary 2-club.)");
			System.out.println("           \t (A t-hereditary 2-club with at least t+2 vertices is a (t+1)-connected 2-club.)");
			System.out.println("-out PREFIX\t Write solution(s) to output file(s) starting with string PREFIX.");
			System.out.println("           \t Program will append value of parameter t and ending '.dimacs'.");
			System.out.println("           \t For each choice of the parameter t one output file will be generated.");
			System.out.println("           \t Output files are in DIMACS format and located in directory 'output'.");
			System.out.println("           \t Default: No output files will be generated.");
			System.out.println("-t X       \t Run algorithm for parameter t = X.");
			System.out.println("           \t Note: If you want to compute a X-hereditary 2-club you have to type '-t X+1'.");
			System.out.println("           \t (See relationship of t-robust and (t-1)-hereditary 2-clubs. ");
			System.out.println("-t MIN-MAX \t Run algorithm for parameter values MIN to MAX.");
			System.out.println("           \t Default: Run algorithm for every possible value of t.");
			System.out.println("-min X     \t Only look for solutions of size at least X.");
			System.out.println("-metis     \t Tell program that input is a Metis file.");
			System.out.println("-dimacs    \t Tell program that input is a DIMACS file.");
			System.out.println("-clq       \t Tell program that input is a DIMACS clq file.");
			System.out.println("           \t The program normally guesses the correct file type by its ending.");
			System.out.println("-silent    \t Generate no console output whatsoever when running algorithm.");
			System.out.println("           \t (only makes real sense when generating output files)");
			System.out.println("-verbose   \t Generate more detailed console output when running algorithm.");
			System.out.println("-nokernels \t Tells the program to NOT kernelize the graph, but try to solve it all at once.");
			System.out.println("-drop X    \t Tells the program to drop the first X kernels during each search.");
			System.out.println("           \t Useful when the algorithm was previously run for the first X kernels.");
			System.out.println("-timeout X \t Sets time limit for each value of t to X seconds.");
			System.out.println("           \t Other algorithms, as for instance those used for graph analysis");
			System.out.println("           \t are also subject to this timeout. Timed out algorithms are");
			System.out.println("           \t not guaranteed to produce optimal solutions.");
			System.out.println("-aX   	   \t Analysis level of the algorithm. How much information about the input graph");
			System.out.println("           \t and solutions will be generated. Replace X by 0-4. Default is 1.");
			System.out.println("           \t Level 0: Do not analyze the graphs.");
			System.out.println("           \t Level 1: Only generate basic information computable in linear time.");
			System.out.println("           \t Level 2: Generate additional information requiring higher order polynomial time.");
			System.out.println("           \t Level 3: Most information. Warning: Requires solving NP-hard problems.");
			System.out.println("           \t Results of the analysis can be found in the headers of generated output files.");
			System.out.println("           \t Furthermore analysis results of the input will be output to the console if not");
			System.out.println("           \t using -silent, as well as results for the solutions if using -verbose.");
			System.out.println();
			System.out.println("To run a correctness self test of the algorithm:   java RobustTwoClub -selftest");
			System.out.println("To run a self test for t-robust variant:           java RobustTwoClub -selftestvb");
			System.out.println("To run a self test for t-connected variant:        java RobustTwoClub -selftesttc");
			System.out.println("To show this usage information:                    java RobustTwoClub -help\n");
			System.exit(0);
		}
		
		if (args[0].equals("-selftest")) {
			System.out.println("Running algorithm selftest for correctness.");
			CorrectnessChecker.runErdosRenyiTests(17000, 4, 20, 0.1, 0.4, 1, 20, true, Model.HEREDITARY);
			System.exit(0);
		} else if (args[0].equals("-selftestvb")) {
			System.out.println("Running algorithm selftest for correctness. (using Veremyev/Boginski model)");
			CorrectnessChecker.runErdosRenyiTests(17000, 4, 20, 0.1, 0.4, 1, 1, true, Model.VB_MODEL);
			System.exit(0);
		} else if (args[0].equals("-selftesttc")) {
			System.out.println("Running algorithm selftest for correctness. (using biconnected model)");
			CorrectnessChecker.runErdosRenyiTests(17000, 4, 20, 0.1, 0.4, 1, 20, true, Model.BICONNECTED);
			System.exit(0);
		} else if (args[0].equals("-selftestall")) {
			CorrectnessChecker.clusteringInstanceTest();
			for (int i = 4; i < 25; i++) {
				System.out.println("Running algorithm selftest for correctness.");
				int tests = 10000;
				CorrectnessChecker.runErdosRenyiTests(tests, i, i, 0.1, 0.4, 1, 20, true, Model.HEREDITARY);
				System.out.println("Running algorithm selftest for correctness. (using Veremyev/Boginski model)");
				CorrectnessChecker.runErdosRenyiTests(tests, i, i, 0.1, 0.4, 1, 20, true, Model.VB_MODEL);
				System.out.println("Running algorithm selftest for correctness. (using biconnected model)");
				CorrectnessChecker.runErdosRenyiTests(tests, i, i, 0.1, 0.4, 1, 20, true, Model.BICONNECTED);
			}
			System.exit(0);
		} 	
		
		String filename = args[0];
		boolean typeGiven = false;
		boolean tSpecified = false;
		boolean noKernels = false;
		Model usedModel = Model.HEREDITARY;
		boolean tb = false;
		int type = 0;					// dummy initialization
		int minT = 1;					// dummy initialization
		int maxT = 1;					// dummy initialization
		int minSize = 0;
		int drop = 0;
		int timeout = -1;
		int verboseness = 1;
		int analysis = 1;
		boolean output = false;
		String outputPrefix = null;
		
		if (args.length > 1) {
			boolean wellFormed = true;
			for (int arg = 1; arg < args.length; arg++) {
				try {
					switch (args[arg]) {
						case "-silent" : 	verboseness = 0; break;
						case "-verbose" : 	verboseness = 2; break;
						case "-a0"	:		analysis = 0; break;
						case "-a1"	:		analysis = 1; break;
						case "-a2"	:		analysis = 2; break;
						case "-a3"	:		analysis = 3; break;
						case "-nokernels" :	noKernels = true; break;
						case "-dimacs" : 	type = 0; typeGiven = true; break;
						case "-metis" : 	type = 1; typeGiven = true; break;
						case "-clq" : 		type = 2; typeGiven = true; break;
						case "-vb" : 		usedModel = Model.VB_MODEL; break;
						case "-tb" : 		tb = true; break;
						case "-tc" : 		usedModel = Model.BICONNECTED; break;
						case "-out" : 		output = true;
											if (args.length < arg) wellFormed = false;
											else outputPrefix = args[++arg];
											break;
						case "-min" :		if (args.length < arg) wellFormed = false;
											else minSize = Integer.parseInt(args[++arg]);
											break;
						case "-timeout" :	if (args.length < arg) wellFormed = false;
											else timeout = Integer.parseInt(args[++arg]);
											break;
						case "-drop" :		if (args.length < arg) wellFormed = false;
											else drop = Integer.parseInt(args[++arg]);
											break;
						case "-t" :			if (args.length < arg) wellFormed = false;
											else {
												tSpecified = true;
												String[] values = args[++arg].split("-");
												if (values.length > 2) { wellFormed = false; break; }
												if (values.length == 1) {
													int t = Integer.parseInt(values[0]);
													if (t < 1) {
														System.out.println("Specified parameter t is invalid.");
														return;
													}
													minT = t; maxT = t;
												} else {
													minT = Integer.parseInt(values[0]);
													maxT = Integer.parseInt(values[1]);
													if (minT > maxT || minT < 1) {
														System.out.println("Specified range for parameter t is invalid.");
														return;
													}
												}
											}
											break;
						default :			wellFormed = false;
					} 
				} catch (NumberFormatException e) {
					wellFormed = false; break;
				}
				if (!wellFormed) break;
			}
			if (!wellFormed) {
				System.out.println("Program does not understand arguments of program call.");
				System.out.println("Please call 'java RobustTwoClub -help' for usage information.");
				return;
			}
		}	

		System.out.println();
		
		if (!typeGiven) type = RtcGraph.guessInputType(filename);
		RtcGraph graph = new RtcGraph(filename, type);
		if (!tSpecified) {
			minT = 1;
			maxT = graph.size();
		}
		RobustTwoClubAlgorithm.run(graph, filename, minT, maxT, minSize, verboseness, analysis, output, outputPrefix, 
					   noKernels, drop, timeout, tb, usedModel);
		long time = System.currentTimeMillis() - start;
		System.out.println("\nTotal program running time: "+(time/1000.0)+" seconds");
	}
}	
