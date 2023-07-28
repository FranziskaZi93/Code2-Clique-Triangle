package robustTwoClub.algorithms;


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.ListIterator;
import java.util.Stack;


import robustTwoClub.analysis.RtcGraphAnalysis;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.latex.TableMaker;

public class RobustTwoClubAlgorithm {
	
	public enum Model {
	    HEREDITARY, 	// t-hereditary 2-Club model, where each pair of non-adjacent vertices needs t common neighbors
	    VB_MODEL, 		// Robust 2-Club model of Veremyev and Boginski, where adjacent vertices also need t common neighbors.
	    BICONNECTED		// biconnected 2-Club model, where the solution has to biconnected and a 2-Club
	}
	
	private static int verbosenessLevel;  // verboseness level of algorithm, i.e. how much information
                                              // is output to console during operation, see description of run()
	private static int analyzeLevel;      // analysis level of algorithm, i.e. how much information about the
					      // input graph and solutions shall be generated, see description of run()
	private static HashSet<Integer> clique;			// Vertices of largest found clique by analysis of input graph
	private static boolean maxClique;				// Whether largest found clique is guaranteed to be maximal
	
	private static String graphName;				// Name of the input graph, usually name of the input file.
	private static String[] headerInfo;				// String array for storing header information for output files.
	
	private static boolean output;					// Whether algorithm should write solutions to files
	private static String outputPrefix;				// Name prefix for generated output files.
	private static int minSolutionSize;
	private static int dropKernels;
	private static Model usedModel;					// Remembers which model is used when runing of the algorithm 

	private static boolean useTriangleBound;
	private static int timeout;
	private static long start;
	private static boolean timeRunOut;
	
	private static int initialSize;					// Initial size of original input graph
	
	private static int t;							// parameter t of problem, i.e. required # of common neighbors
	private static int lastBestT;					// Remembers maximum t for which last solution is a 2,t-club
													// to skip values of t for which the last solution applies.
	private static RtcGraph solution;				// stores currently best solution 
	private static int lowerBound;					// size of currently best solution as lower bound
	
	private static ArrayList<Integer> kernelList;	// List of kernels sorted by size of 2-neighborhood in ascending order
	private static int[] kernelSizes;
	
	private static RtcGraph kernel;					// the current state of the currently considered kernel
	private static HashSet<Integer> marked;			// the currently marked vertices of the kernel
	private static int[][] noCommonNeighbors;		// # common neighbors for all vertex pairs of kernel,
													// triangular matrix only containing pair (w,v) for IDs v < w, not (v,w)
	private static boolean[][] agree;				// true if the two vertices satisfy conditions (being at distance at most 2 plus side-constraint) of respective model in the graph 'kernel' 
	private static int[][] disagreeReason;			// if two vertices u,v where agreeable once then disagreeReason contains the vertex whose deletion made u,v disagreeing 
													// (when deleting disagreeReason[u][v] the value of agree[u][v] was set from true to false)
													// only used for model BICONNECTED
	private static int[][] connectednessLB;			// a lower bound for the number of internally vertex disjoint paths between the two vertices 
	private static boolean[][] connectednessReset;	// whether the last value in connectednessLB was computed exactly with the flow algorithm
	private static int[] agreements;				// remembers for every vertex the # of agreements
	
	// remembers vertex pairs in conflict with each other (to speed up VCR), stores each conflict (v, w) once for v and w
	private static HashMap<Integer, HashSet<Integer>> conflicts;
	
	private static Stack<Integer> rollbackStack;	// maintains order of vertex operations to ensure LIFO order undo
													
	private static final DecimalFormat tf = new DecimalFormat("0.000");		// Format for time output
	private static final DecimalFormat df = new DecimalFormat("00.00");		// Format for progress output
	
	private static int depth, branches, ops, completedOps, hardKernels, initializedKernels, largestKernel, 
					largestHardKernel, vcrApplications,	kernelSizeSum, hardKernelSizeSum;
	private static long	vcrTime, branchTime, opsTime, kernelBuildTime, timeSpentOnMaxFlow,
					graphReductionTime, kernelReductionTime, kernelInitializationTime, initialReductionTime;
	private static long maxBranches;
	private static long curBranches;
	private static long maxFlowCalls, posAnswer, negAnswer, conflictChecks;
	private static long[] maxFlowCallsPerStage;
	private static int stage;
	
	
	
	private static double R2KernelReductionSum; // sum of reduction percentages of initial R2/R3 on kernels
	private static double R2KernelReductionMin; // minimum reduction achieved by R2/R3
	private static double R2KernelReductionMax; // maximum reduction achieved by R2/R3
	private static int R2SurvivingKernels;		// total number of kernels not dismissed WHILE performing R2/R3
	private static double R4KernelReductionSum; // sum of reduction percentages of initial R4 on kernels
	private static double R4KernelReductionMin; // minimum reduction achieved by R4
	private static double R4KernelReductionMax; // maximum reduction achieved by R4
	private static int R4SurvivingKernels;		// total number of kernels not dismissed WHILE performing R4
	
	// The following information allows us to eliminate superfluous applications of the Vertex Cover Rule
	private static int newConflicts;				// Counts new conflicts since the last application of the VCR,
													// new conflicts covered by the last VC are not counted
	private static HashSet<Integer> lastVC;			// Remembers the last vertex cover of the conflict graph,
													// so we only count new conflicts not covered by the old VC
	
	public static long solveTime;					// Time needed (in msecs) by the last call of the algorithm
	public static int solutionSize;					// Size of the last solution generated by the algorithm
	public static int neededBranches;				// Search tree branches generated by the last run of the algorithm
	public static int neededOps;					// Operations required by the last run of the algorithm
	public static boolean timeoutOnLastT;
	public static int lastT;
	
	// Dirty solution for all the old method calls without the 'vb' parameter. 
	public static RtcGraph run (RtcGraph graph, String name, int minT, int maxT, int minSize, int vLevel, 
				    int aLevel,	boolean out, String prefix, boolean noKernels, int drop, int timeOut,boolean tb) {
	    return run(graph, name, minT, maxT, minSize, vLevel, aLevel, out, prefix, noKernels, drop, timeOut, tb, Model.HEREDITARY);
	}
	
	/**
	 * 
	 * @param graph		The graph on which 2,t-clubs shall be found. Graph will be reduced by algorithm!
	 * 					So pass a copy, if you still need the original graph!
	 * @param name		Name of the input graph, usually name of input file. Pass empty string if name not required, not null!
	 * @param minT		Smallest value of parameter t to check.
	 * @param maxT		Largest value of parameter t to check, algorithm aborts when it can't find solution for some value,
	 * 					because then larger values of t are guaranteed to produce no (non-trivial) solution.
	 * @param minSize	Required minimum size of solutions. 
	 * @param vLevel	Verboseness level.		0 = Be quiet! No console output whatsoever. 
	 *											1 = Just basic information and progress report.
	 *											2 = Detailed informations about run of algorithm.
	 *											3 = Detailed informations, but no progress report.
	 * @param aLevel	Analysis level, how much the input graph and solutions shall be analyzed.
	 * 											0 = Do not analyze the input graph and solutions.
	 *											1 = Only basic information computable in linear time.
	 * 											2 = Additional information requiring higher order polynomial time.
	 * 											3 = Even more information requiring to solve NP-hard problems.
	 * @param out		Whether algorithm shall write solutions to disk.
	 * @param prefix	Prefix for name output files. Just pass null if 'output' is set to false.
	 * @param noKernels	Whether kernelization should be deactivated. Kernelization is necessary for larger graphs
	 * 					(more than a few thousand vertices or more than about ten to twenty thousand edges)!
	 * @param drop		Drop X first kernels during search.
	 * @param timeOut	Allowed time in seconds for each run of a parameter t. Set to -1 for no limit.
	 * @param vb		Whether the robust 2-club model of Veremyev/Boginski shall be used instead,
	 * 					where two adjacent vertices also require t common neighbors.
	 * @return			Solution for parameter maxT, null if none.
	 * 					(Returning solution is required by correctness checker.)
	 * @throws			GraphException
	 * @throws 			AlgorithmException
	 */
	public static RtcGraph run (RtcGraph graph, String name, int minT, int maxT, int minSize, int vLevel, int aLevel,
				    boolean out, String prefix, boolean noKernels, int drop, int timeOut, boolean tb, Model vb) {
		RtcGraph copy = graph.getClone();
		graphName = name;
		verbosenessLevel = vLevel;
		analyzeLevel = aLevel;
		output = out;
		outputPrefix = prefix;
		minSolutionSize = minSize;
		dropKernels = drop;
		useTriangleBound = tb;
		usedModel = vb;
		timeout = timeOut;
		solution = null;			// Purge solution from previous call of algorithm
		solutionSize = 0;
		lastBestT = 0;
		initialSize = graph.size();
		
		timeSpentOnMaxFlow = 0;
		maxFlowCalls = 0; posAnswer = 0; negAnswer = 0; conflictChecks = 0;
		maxFlowCallsPerStage = new long[10];
		stage = 0;
		RtcGraph.flowComputePath = 0;
		RtcGraph.flowConstructGraph = 0;
		
		// Generate some basic information on graph, if analyze level isn't set to 0.
		clique = null;
		if (analyzeLevel == 3) clique = new HashSet<Integer>();
		maxClique = false;
		if (analyzeLevel > 0 && (verbosenessLevel > 0 || output)) {
			maxClique = RtcGraphAnalysis.analyzeGraph(graph, analyzeLevel, graphName, verbosenessLevel < 2, timeout, false);
			clique = RtcGraphAnalysis.clique;
			headerInfo = RtcGraphAnalysis.header;
			if (verbosenessLevel > 1) System.out.println();
		} else headerInfo = null;
		
		// Initialize LaTeX TableMaker, if output = true and analyze level = 3
		if (output && analyzeLevel == 3) TableMaker.init(graphName);
		
		long start = System.currentTimeMillis();
		
		// Build kernel list by sorting vertices of graph by size of 2-neighborhood in ascending order
		kernelSizes = new int[graph.size()];

		initialReductionTime = 0;
		if (minT > 1)
		{
			t = minT;
			long time = System.nanoTime();
			// Shrink original graph by applying R2 (deleting low degree vertices)
			HashSet<Integer> deleteSet = new HashSet<Integer>();
			int tmpSize = graph.size();
			for (int vertex : graph.getVertices())
			{
				if (graph.degree(vertex) < t)
					deleteSet.add(vertex);
			}
			for (int vertex : deleteSet)
				if (graph.contains(vertex))					// recursive deletions might have already deleted it
					deleteGraphVertex(graph, vertex, null, true);
			initialReductionTime += System.nanoTime() - time;
			if (verbosenessLevel > 0) System.out.println("Preprocessing: Removed "+ (tmpSize - graph.size()) +" low-degree vertices. ");
			t = 0;
		}
		
		RtcGraph returnSolution = null;
		RtcGraphAnalysis.setCompareSolution(null);
		
		if (graph.size() > 0)
		{
			// Build kernel list by sorting vertices of graph by size of 2-neighborhood in ascending order
			if (!noKernels) {
				if (verbosenessLevel > 0) System.out.println("Preprocessing: Building kernel list.");
				ArrayList<Integer> unsorted = new ArrayList<Integer>();
				for (int id : graph.getVertices()) {
					kernelSizes[id] = graph.sizeOfTwoNeighborhood(id, true);
					unsorted.add(id);
				}
				kernelList = kernelQuicksort(unsorted);
				if (verbosenessLevel > 0)
					System.out.println("Required preprocessing time: " + 
							tf.format((System.currentTimeMillis() - start) / 1000.0) + " seconds.");
			} else for (int v = 0; v < graph.size(); v++)
				kernelSizes[v] = graph.size();	// Dummy values to avoid deletions because of "kernel size"
			
			// Run algorithm for all choices of the parameter t specified by range [minT; maxT]
			for (t = minT; t <= maxT; t++) {
				// Check whether last solution is also a solution for higher t, so we can skip some steps
				if (lastBestT >= t) {
					if (verbosenessLevel > 0) 
						System.out.println("\nLargest 2,"+(t-1)+"-club is also a 2,"+lastBestT+"-club.");
					t = lastBestT;
					continue;
				}
				/* Run algorithm for current parameter t. Stop algorithm as soon as any value of the parameter
				 * t yields no solution, greater values of t are then guaranteed to yield no solution.
				 */
				returnSolution = runForT(graph, noKernels, (maxT > t));
				if (returnSolution == null && !timeRunOut) break;
				// Check whether found solution is a clique.
				if (!timeRunOut && returnSolution.isClique()) {
					if (verbosenessLevel > 0) {
						System.out.println("\nLargest 2,"+t+"-club is a clique. "+
								"This is also the best solution for all t > "+t+".");
						System.out.println("This size "+returnSolution.size()+
								" clique is also the largest clique of the graph.");
					} break;
				}
				RtcGraphAnalysis.setCompareSolution(solution);
			}
		}
		if (verbosenessLevel > 0) System.out.println();
		// If last run returned no solution we search for the maximum clique (which is the best solution)
		// (in the Veremyev/Boginski or BICONNECTED model small cliques are NOT solutions)
		if (returnSolution == null && minSize <= t && !timeRunOut && usedModel == Model.HEREDITARY) { 
			if (clique == null) {
				if (verbosenessLevel > 0) 
					System.out.println("Running clique algorithm to determine solution for t > " + (t - 1));
				clique = CliqueAlgorithm.run(copy.getClone(), timeout, verbosenessLevel == 0);
				maxClique = !CliqueAlgorithm.timedOut;
			} else if (verbosenessLevel > 0)
				System.out.println("Using clique determined by initial analysis as solution for t > " + (t - 1));
			if (clique == null)
				System.out.println("Clique algorithm timed out without producing any solution.");
			else { 
				if (!maxClique && verbosenessLevel > 0)
					System.out.println("The clique algorithm timed out, so this is clique is not guaranteed to be maximum.");
				solution = copy;
				solution.retainVerticesByID(clique);
				if (verbosenessLevel > 0)
					System.out.println("Maximum clique has size "+solution.size()+".");
				String[] solutionHeader = new String[1];
				if (maxClique) solutionHeader[0] = "This solution is a maximum clique of the original graph and has size "+
						solution.size()+".";
				else solutionHeader[0] = "This solution is a clique and has size " + solution.size() + ".";
				String[] header = RtcGraphAnalysis.combineHeaders(headerInfo, solutionHeader, t, graphName, maxClique, timeout);
				if (output) solution.writeToFile("output", outputPrefix+t, header);
				if (output && verbosenessLevel > 0)
					System.out.println("Found clique written to '"+outputPrefix+t+".dimacs' in output folder.");
				// Info for TableMaker
				if (output && analyzeLevel == 3) TableMaker.addCliqueLine(t, solution.size());
				returnSolution = solution.getClone();
			}
		}
		
		if (output && analyzeLevel == 3) TableMaker.finalizeAndOutput("output", outputPrefix+"table");
		
		solveTime = System.currentTimeMillis() - start;
		solutionSize = 0;
		neededBranches = branches;
		neededOps = ops;
		lastT = t;
		if (solution != null) solutionSize = solution.size();
		return returnSolution;	// Return solution for last value of parameter t. (correctness checker needs this)
	}	
	
	private static ArrayList<Integer> kernelQuicksort (ArrayList<Integer> toSort) {
		int pivotElem = toSort.get((int) (Math.random() * (toSort.size())));
		int pivotVal = kernelSizes[pivotElem];
		ArrayList<Integer> equalsList = new ArrayList<Integer>();
		ArrayList<Integer> smallerList = new ArrayList<Integer>();
		ArrayList<Integer> greaterList = new ArrayList<Integer>();
		for (int elem : toSort) {
			if (kernelSizes[elem] < pivotVal) smallerList.add(elem);
			else if (kernelSizes[elem] > pivotVal) greaterList.add(elem);
			else equalsList.add(elem);
		}
		ArrayList<Integer> solution = new ArrayList<Integer>();
		if (!smallerList.isEmpty()) solution.addAll(kernelQuicksort(smallerList));
		solution.addAll(equalsList);
		if (!greaterList.isEmpty()) solution.addAll(kernelQuicksort(greaterList));
		return solution;
	}
	
	private static RtcGraph runForT (RtcGraph g, boolean noKernels, boolean needPrivateCopy) {
		
		start = System.currentTimeMillis();
		timeRunOut = false;
		
		// Reset lower bound and solution
		lowerBound = 0;
		solution = null;
		
		// Reset reduction rule statistics
		R2KernelReductionSum = 0.0;
		R2KernelReductionMin = 100.0;
		R2KernelReductionMax = 0.0;
		R2SurvivingKernels = 0;
		R4KernelReductionSum = 0.0;
		R4KernelReductionMin = 100.0;
		R4KernelReductionMax = 0.0;
		R4SurvivingKernels = 0;
		
		if (verbosenessLevel > 0) {
			System.out.println("\nRunning search for largest 2,"+t+"-club...");
		}
		
		if (t == 1) {// use closed neighborhood of max degree vertex as initial solution for t = 1
			int maxDegreeVertex = -1;
			for (int v : g.getVertices())
				if (maxDegreeVertex == -1 || g.degree(v) > g.degree(maxDegreeVertex))
					maxDegreeVertex = v;
			if (! (usedModel == Model.BICONNECTED && g.degree(maxDegreeVertex) == 0) )
			{
				solution = g.getOneNeighborhood(maxDegreeVertex);
				lowerBound = solution.size();
			}
		}
		
		if (clique != null && lowerBound < clique.size() && usedModel == Model.HEREDITARY) {
//			solution = null;
			solution = g.getSubgraph(clique, clique.iterator().next());
			lowerBound = clique.size();
		}	
		
		// Set lower bound to required solution size and throw away approximate initial solution if too small
		if (lowerBound < minSolutionSize) {
			lowerBound = minSolutionSize - 1; 
			solution = null;
			solutionSize = 0;
		}
		
		if (verbosenessLevel > 1 && lowerBound > 0)
			System.out.println("Initial lower bound: "+ lowerBound);
		
		if (verbosenessLevel > 0 && verbosenessLevel != 3)
			if (!noKernels) System.out.print("Preprocessing kernels...");
		
		graphReductionTime = initialReductionTime;
		long time = System.nanoTime();
		// Shrink original graph by applying R2 (deleting low degree vertices)
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int vertex : g.getVertices())
			if (g.degree(vertex) < t || g.degree(vertex) == 1)
				deleteSet.add(vertex);
		for (int vertex : deleteSet)
			if (g.contains(vertex))					// recursive deletions might have already deleted it
				deleteGraphVertex(g, vertex, kernelSizes, true);
		graphReductionTime += System.nanoTime() - time;
//		System.out.println("");
//		System.out.println("Graph reduction time 1: " + tf.format(graphReductionTime/1000000000.0) + " seconds" + "\t\t #vertices: " + deleteSet.size());

		// Remember reduction percentage of initial application of R2 on graph
		int reducedSize = g.size();
		double reductionPercentage = (initialSize - reducedSize) * 100.0 / initialSize;
		
		// After this point we will start making changes to the graph which will invalidate it for subsequent
		// runs of the algorithm for other values of t. Therefore we create a work copy of it, if we make further runs.
		RtcGraph graph = g;
		int[] kernelSizesLocal = kernelSizes;
		if (needPrivateCopy) {
			graph = g.getClone();
			kernelSizesLocal = kernelSizes.clone();
		}
		
		// For t=1 (when not using the VB model) or a specified minimum solution size, we have an 
		// initial lower bound and can remove kernels whose size is below it.
		time = System.nanoTime();
		if (lowerBound > 0 && !noKernels) {
			deleteSet = new HashSet<Integer>();
			for (int vertex : graph.getVertices())
				if (kernelSizesLocal[vertex] <= lowerBound)
					deleteSet.add(vertex);
			for (int vertex : deleteSet)
				if (graph.contains(vertex))						// recursive deletions might have already deleted it
					deleteGraphVertex(graph, vertex, kernelSizesLocal, false);	
		}
		graphReductionTime += System.nanoTime() - time;

		int kernelsSolved = 0;
		if (!noKernels && verbosenessLevel > 0 && verbosenessLevel != 3) {
			System.out.print("\rProgress: 00.00%    ");
			System.out.print("Current Lower Bound: " + lowerBound + " \t");
			System.out.print("Vertices Remaining: " + graph.size()+ " \t");
			System.out.print("Kernels: " + kernelList.size());
		}
		depth = 0; branches = 0; ops = 0; completedOps = 0; hardKernels = 0; largestKernel = 0; largestHardKernel = 0;
		vcrApplications = 0; vcrTime = 0; branchTime = 0; opsTime = 0; kernelBuildTime = 0; kernelReductionTime = 0;
		kernelInitializationTime = 0; kernelSizeSum = 0; hardKernelSizeSum = 0; initializedKernels = 0;
		maxBranches = 0;
		long timer = System.currentTimeMillis();		// Used for generating progress report every second
		// If kernelization is active, consider kernels in order of ascending original size
		int dropCounter = dropKernels;
		if (!noKernels) {
			int kernels = kernelList.size();
			int kernelsChecked = 0;
			ListIterator<Integer> iterator = kernelList.listIterator();
			while (iterator.hasNext()) {
				if (timeout != -1 && System.currentTimeMillis() - start > timeout*1000) {
					timeRunOut = true; break;
				}
				int center = iterator.next();
				if (dropCounter <= 0 && graph.contains(center) && kernelSizesLocal[center] > lowerBound) {
					long timeStamp = System.nanoTime();
					// kernelSizesLocal contains only upper bounds after deletions, so check actual kernel size!
					int actualSize = graph.sizeOfTwoNeighborhood(center, true);
					if (actualSize  <= lowerBound) continue;
					if (actualSize > largestKernel)
						largestKernel = actualSize;
					kernel = graph.getTwoNeighborhood(center);
					kernelBuildTime += System.nanoTime() - timeStamp;
					kernelSizeSum += kernel.size();
					searchOnKernel(true);		// Search on kernel with kernel center vertex forced to being marked
					stage = 0;
					kernelsSolved++;
				}
				time = System.nanoTime();
				if (graph.contains(center)) deleteGraphVertex(graph, center, kernelSizesLocal, false);
				graphReductionTime += System.nanoTime() - time;
				kernelsChecked++;
				dropCounter--;
				if (verbosenessLevel > 0 && verbosenessLevel != 3 && System.currentTimeMillis() - timer > 1000) {
					timer = System.currentTimeMillis();
					System.out.print("\rProgress: " + df.format((double)(100*kernelsChecked)/kernels) + "%    ");
					System.out.print("Current Lower Bound: " + lowerBound + " \t");
					System.out.print("Vertices Remaining: " + graph.size() + " \t");
					if (verbosenessLevel > 1 && usedModel == Model.BICONNECTED)
					{
						System.out.print("Max-Flow: "+df.format(timeSpentOnMaxFlow/1000000000.0) + "s (" + df.format((System.currentTimeMillis() - start)/1000.0) + "s total) \t [");
						System.out.print("Conflic checks: " +conflictChecks + "  Flow calls: "+ maxFlowCalls + " pos/neg:" + posAnswer + "/" + negAnswer +" " +Arrays.toString(maxFlowCallsPerStage) + " \t");
						System.out.print("Init Graph: "+df.format(RtcGraph.flowConstructGraph/1000000000.0) + "s \t");
						System.out.print("Ford-Fulkerson: "+df.format(RtcGraph.flowComputePath/1000000000.0) + "s] \t");
						System.out.print("ops: " + ops + " \t");
						System.out.print("branches: " + branches + " \t");
					}
				}
			}
		} else {	// If noKernels is set try to solve the whole graph in one run
			long timeStamp = System.nanoTime();
			kernel = graph.turnIntoSingleKernel();
			kernelBuildTime += System.nanoTime() - timeStamp;
			largestKernel = graph.size();
			searchOnKernel(false);		// Search on single kernel without marking any vertex initially
			kernelsSolved = 1;
		}
		double solveTime = (System.currentTimeMillis() - start);
		if (verbosenessLevel > 0 && timeRunOut) {
			System.out.print("\rTimeout of "+timeout+" seconds was reached.");
			System.out.println("                                                                   ");
		} else if (verbosenessLevel > 0) {
			System.out.print("\rTime required: "+tf.format(solveTime / 1000.0) + " seconds.");
			System.out.println("                                                                   ");
		}
		if (verbosenessLevel > 1) {
			System.out.println("Kernels: "+kernelsSolved+"  Initialized kernels: "+initializedKernels+
					"  Hard kernels: "+hardKernels);
			System.out.println("Largest kernel: "+largestKernel+"  Largest hard kernel: "+largestHardKernel+
					"  Avg. kernel size: "+kernelSizeSum/Math.max(1, kernelsSolved)+
					"  Avg. hard kernel size: "+hardKernelSizeSum/Math.max(1, hardKernels));
			System.out.println("Max. branch depth: "+depth+"  Total branches: "+branches+"  Average branches: " +
					((hardKernels>0)?(branches/hardKernels):0) + "  Max. branches: " + maxBranches + 
					"  Ops: "+ ops + "  Completed Ops: "+completedOps+"  VCR applications: "+vcrApplications);
			System.out.println("Branch time: " + tf.format(branchTime/1000000000.0) + " seconds" +
					"  Ops time: " + tf.format(opsTime/1000000000.0) + " seconds" +
					"  VCR time: " + tf.format(vcrTime/1000000000.0) + " seconds" +
					"  Kernel initialization time: " + tf.format(kernelInitializationTime/1000000000.0) + " seconds");
			System.out.println("Graph reduction time: " + tf.format(graphReductionTime/1000000000.0) + " seconds"+
					"  Kernel build time: " + tf.format(kernelBuildTime/1000000000.0) + " seconds"+
					"  Kernel reduction time: " + tf.format(kernelReductionTime/1000000000.0) + " seconds");
			System.out.println("Conflict checks: " + conflictChecks);
			if (usedModel == Model.BICONNECTED)
			{
				System.out.println("\nTime spent on Max-Flow-Computations: "+df.format(timeSpentOnMaxFlow/1000000000.0) + " seconds from " + df.format((System.currentTimeMillis() - start)/1000.0) + " seconds overall so far.");
				System.out.println("Flow time contains (for "+maxFlowCalls+" calls:");
				System.out.println("Graph construction time: "+df.format(RtcGraph.flowConstructGraph/1000000000.0) + " seconds.");
				System.out.println("Path finding time: "+df.format(RtcGraph.flowComputePath/1000000000.0) + " seconds.");
				System.out.println();
			}
		}	
		
		// Report solution and output it, if output enabled
		if (solution != null && verbosenessLevel > 0) {
			System.out.println("Found solution of size " + lowerBound + " with " +
					solution.getEdgeCount() + " edges.");
			if (timeRunOut)
				System.out.println("Solution is NOT GUARANTEED TO BE OPTIMAL due to timeout.");
			else if (dropKernels > 0)
				System.out.println("Note that solution is NOT GUARANTEED TO BE OPTIMAL due to dropped kernels.");
		} else if (verbosenessLevel > 0 && minSolutionSize <= t) {
			if (!timeRunOut && dropKernels == 0) 
				System.out.println("Found no solution. Solution must be clique of size < " +(t+1)+". "+
					"This also applies for all t > "+t+".\n"+
					"Graph is guaranteed to contain no clique of size > "+t+".\n");
			else if (timeRunOut) System.out.println("Found no solution within allowed time limit.\n");
			else System.out.println("Found no solution in kernels searched.\n");
		} else if (verbosenessLevel > 0) {
			if (!timeRunOut && dropKernels == 0) 
				System.out.println("Found no solution of size at least "+minSolutionSize+".\n");
			else if (timeRunOut)
				System.out.println("Found no solution of size at least "+minSolutionSize+" within allowed time limit.\n");
			else System.out.println("Found no solution of size at least "+minSolutionSize+" in kernels searched.\n");
		}
			
		if (solution != null && output) {
			// Generate header information
			String[] solutionHeader = null;
			// analyzeSolution will also check whether our solution is a solution for higher t, so that we can skip steps
			if (solution.isClique()) {
				solutionHeader = new String[1];
				if (!timeRunOut) solutionHeader[0] = "This solution is a maximum clique of the original graph and has size "+
						solution.size()+".";
				else solutionHeader[0] = "This solution is a clique and has size " + solution.size() + ".";
				lastBestT = solution.is2TClub(usedModel);
			} else if (analyzeLevel > 0) {
				if (verbosenessLevel > 1) System.out.println();
				lastBestT = RtcGraphAnalysis.analyzeSolution(solution, analyzeLevel, t, verbosenessLevel < 2, 
						timeout, output, usedModel);
				solutionHeader = RtcGraphAnalysis.header;
				// Pass information to TableMaker
				if (analyzeLevel == 3) {
					TableMaker.addAlgorithmInfo(solveTime/1000.0, branchTime/1000000000.0, kernelsSolved,
							hardKernels, largestHardKernel, branches, ops);
				}
				if (verbosenessLevel > 1) System.out.println();
			} else lastBestT = solution.is2TClub(usedModel);
			// Combine input header and solution header to single header
			String[] header = RtcGraphAnalysis.combineHeaders(headerInfo, solutionHeader, t, graphName, !timeRunOut, timeout);
			// Write solution to file
			solution.writeToFile("output", outputPrefix+t, header);
			if (verbosenessLevel > 0)
				System.out.println("Found 2,"+t+"-club written to '"+outputPrefix+t+".dimacs' in output folder.");
		} else if (solution != null) {
			if (verbosenessLevel > 1 && analyzeLevel > 0) {
				System.out.println();
				lastBestT = RtcGraphAnalysis.analyzeSolution(solution, analyzeLevel, t, false, 
						timeout, output, usedModel);
				System.out.println();
			} else lastBestT = solution.is2TClub(usedModel);
		}
		
		// If verbose flag was used, additionally output information about reduction rule effectiveness
		if (verbosenessLevel > 1) {
			System.out.println("\nAdditional information regarding effectiveness of reduction rules for t="+t+": ");
			System.out.println("Initial application of low degree rule reduced graph from " + initialSize +
				" to " + reducedSize + " vertices (by " + df.format(reductionPercentage) + "%).");
			System.out.println("Out of " + kernelsSolved + " kernels, " + (kernelsSolved - R2SurvivingKernels) + 
					" kernels were dismissed WHILE performing initial pass of R2/R3. " +
					"(" + df.format((kernelsSolved-R2SurvivingKernels)*100.0/kernelsSolved) + "%)");
			if (R2SurvivingKernels > 0) {
				System.out.println("The remaining " + R2SurvivingKernels + " kernels were reduced by between " +
						df.format(R2KernelReductionMin) + "% and " + df.format(R2KernelReductionMax) + "% by R2. "
						+ "(on average " + df.format(R2KernelReductionSum/R2SurvivingKernels) + "%)");
				System.out.println("Out of these " + R2SurvivingKernels + " kernels, " + 
						(R2SurvivingKernels - R4SurvivingKernels) + " kernels were dismissed WHILE performing " +
						"initial marking or the initial application of R4. (" +
						df.format((R2SurvivingKernels-R4SurvivingKernels)*100.0/R2SurvivingKernels) + "%)");
			}
			if (R4SurvivingKernels > 0) {
				System.out.println("The remaining " + R4SurvivingKernels + " kernels were reduced by between " +
						df.format(R4KernelReductionMin) + "% and " + df.format(R4KernelReductionMax) + "% by R4. "
						+ "(on average " + df.format(R4KernelReductionSum/R4SurvivingKernels) + "%)");
				System.out.println("Out of these " + R4SurvivingKernels + " kernels, " + 
						(R4SurvivingKernels - hardKernels) + " kernels were dismissed by the initial " +
						"application of the Vertex Cover Rule (R5). (" +
						df.format((R4SurvivingKernels-hardKernels)*100.0/R4SurvivingKernels) + "%)");
			
			}	
		}
		
		if (timeRunOut) lastBestT = t;	// Reset lastBestT if we had a timeout.
		timeoutOnLastT =  timeRunOut;
		return solution;
	}
	
	// Set onlyByDegree = false to reduce by kernel size and degree, set it to true to reduce by degree only
	private static void deleteGraphVertex (RtcGraph g, int vertex, int sizes[], boolean onlyByDegree) {
		HashSet<Integer> checked = new HashSet<Integer>();
		checked.add(vertex);
		HashSet<Integer> deleteSet = recursiveDelete(g, vertex, sizes, checked, onlyByDegree);
		checked.addAll(deleteSet);
		while (!deleteSet.isEmpty()) {
			int next = 0;								// dummy initialization
			for (int v : deleteSet) {next = v; break;}	// grab any vertex in delete set
			deleteSet.remove(next);
			if (g.contains(next)) {
				HashSet<Integer> newDeletes = recursiveDelete(g, next, sizes, checked, onlyByDegree);
				deleteSet.addAll(newDeletes);
				checked.addAll(newDeletes);
			}	
		}
	}
	
	private static HashSet<Integer> recursiveDelete (RtcGraph g, int vertex, int sizes[], HashSet<Integer> checked, boolean onlyByDegree) {
		HashSet<Integer> seen = new HashSet<Integer>();
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		seen.add(vertex);
		for (int v : g.getNeighbors(vertex))
			if (!checked.contains(v)) {
				if (sizes != null) sizes[v]--;
				seen.add(v);
				if (g.degree(v) <= t || g.degree(v) <= 2) deleteSet.add(v); // R2 will apply after deleting vertex
				if (!onlyByDegree && sizes[v] <= lowerBound) deleteSet.add(v);
			}
		if (sizes != null)
		{	
			for (int v : g.getNeighbors(vertex))
				for (int w : g.getNeighbors(v))
					if (!seen.contains(w) && !checked.contains(w)) {
						sizes[w]--;
						seen.add(w);
						if (!onlyByDegree && sizes[w] <= lowerBound) deleteSet.add(w);
					}
		}
		g.deleteVertex(vertex);
		return deleteSet;
	}
	
	private static void searchOnKernel (boolean markZero) {
		
		int maxID = kernel.size();
		curBranches = 0;
		stage = 1;
		
		long time = System.nanoTime();
		int initialKernelSize = kernel.size();
		// Shrink kernel by applying R2 and R3 (deleting low degree vertices and vertices in conflict with center)
		if (markZero) {		// Single pass variant gets whole graph, which is already reduced w.r.t. vertex degrees
			HashSet<Integer> deleteSet = new HashSet<Integer>();
			for (int vertex : kernel.getVertices()) {
				if (kernel.degree(vertex) < t || kernel.degree(vertex) == 1) {
					if (vertex == 0) return;	// kernel center vertex has too low degree, abort kernel
					deleteSet.add(vertex);		// delete vertex of too low degree
				}
				else if (vertex != 0 && isInConflict(vertex,0,false,0)) 	//!needNoNeighbors(vertex, 0) && kernel.countCommonNeighbors(vertex, 0) < t)
					deleteSet.add(vertex);		// vertex in conflict with kernel center vertex
			}
	        if (kernel.size() - deleteSet.size() <= lowerBound) return;		// Kernel becomes too small, abort
			if (!deleteSet.isEmpty())
				if (deleteKernelVertices(deleteSet)) return;				// Abort if center vertex gets deleted
			if (kernel.size() <= lowerBound) return;						// Kernel has become too small, abort
		}
		kernelReductionTime += System.nanoTime() - time;
		// Collect information for effectiveness of R2/R3 on initial kernels
		int reducedKernelSize = kernel.size();
		double kernelReductionPercentage = (initialKernelSize - reducedKernelSize) * 100.0 / initialKernelSize;
		if (kernelReductionPercentage < R2KernelReductionMin) R2KernelReductionMin = kernelReductionPercentage;
		if (kernelReductionPercentage > R2KernelReductionMax) R2KernelReductionMax = kernelReductionPercentage;
		R2KernelReductionSum += kernelReductionPercentage;
		R2SurvivingKernels++;
		
		stage = 2;

		time = System.nanoTime();
		// Initialize agreeability matrix, agreements vector and conflict set
		noCommonNeighbors = new int[maxID][];
		connectednessLB = new int[maxID][];
		agree = new boolean[maxID][];
		connectednessReset = new boolean[maxID][];
		disagreeReason = new int[maxID][];
		for (int i = 0; i < maxID; i++) 
		{
			noCommonNeighbors[i] = new int[i];
			connectednessLB[i] = new int[i];
			agree[i] = new boolean[i];
			connectednessReset[i] = new boolean[i];
			disagreeReason[i] = new int[i];
		}
		agreements = new int[maxID];
		conflicts = new HashMap<Integer, HashSet<Integer>>();
		for (int v : kernel.getVertices())
			conflicts.put(v, new HashSet<Integer>());
		// Set all entries of agree[v][w] to size of common neighborhood of v and w
		for (int v : kernel.getVertices())
			for (int n1 : kernel.getNeighbors(v))
				for (int n2 : kernel.getNeighbors(v)) if (n1 < n2)
					noCommonNeighbors[n2][n1]++;
		// Increase agreement for agreeable pair and introduce conflicts for non-agreeable pairs
		for (int v : kernel.getVertices())
			for (int w : kernel.getVertices()) if (v < w)
				if  (!isInConflict(w, v, true, 0)) { //(noCommonNeighbors[w][v] >= t || needNoNeighbors(v, w)) {
					agreements[v]++;
					agreements[w]++;
					agree[w][v] = true;
					disagreeReason[w][v] = -1; // init as agreeable
				} else {
					conflicts.get(v).add(w);
					conflicts.get(w).add(v);
					disagreeReason[w][v] = -2; // init as never will be agreeable
				}
		// Initialize rollback stack and marked set
		rollbackStack = new Stack<Integer>();
		marked = new HashSet<Integer>();		
		kernelInitializationTime += System.nanoTime() - time;
		initializedKernels++;
		
		// If kernelization is active mark center vertex of kernel and perform reduction due to initial marking
		time = System.nanoTime();
		if (markZero) if (!mark(0))				// Do not mark vertex 0 if we do not use kernelization
			return;								// mark(0) should actually always return true, but who knows?
		kernelReductionTime += System.nanoTime() - time;
		stage = 3;

		time = System.nanoTime();
		int intermediateKernelSize = kernel.size();
		// Shrink kernel by applying reduction rule R4 (Low Agreement Rule) and checking R2 and R3 again, if R4 applies
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int vertex : kernel.getVertices()) {
			if (agreements[vertex] < lowerBound) {
				if (markZero && vertex == 0)
					return;					// kernel center vertex has too low agreement, abort kernel
				deleteSet.add(vertex);		// delete vertex of too low agreement
			}
		}
		if (kernel.size() - deleteSet.size() <= lowerBound) return;		// Kernel becomes too small, abort
		for (int vertex : deleteSet)
			if (kernel.contains(vertex))								// recursive deletion might have already deleted it
				if (!delete(vertex)) return;
		if (kernel.size() <= lowerBound) return;						// Kernel has become too small, abort
		kernelReductionTime += System.nanoTime() - time;
		
		// Collect information for initial effectiveness of R4/R2/R3 on hard kernels
		reducedKernelSize = kernel.size();
		kernelReductionPercentage = (intermediateKernelSize - reducedKernelSize) * 100.0 / intermediateKernelSize;
		if (kernelReductionPercentage < R4KernelReductionMin) R4KernelReductionMin = kernelReductionPercentage;
		if (kernelReductionPercentage > R4KernelReductionMax) R4KernelReductionMax = kernelReductionPercentage;
		R4KernelReductionSum += kernelReductionPercentage;
		R4SurvivingKernels++;		
		
		stage = 4;
		
		// Initial application of R5 (Vertex Cover Rule)
		int vcrResult = vertexCoverRule(kernel.size() - lowerBound - 1);
		if (kernel.size() - vcrResult <= lowerBound) return;
		newConflicts = 0;
		hardKernels++;
		hardKernelSizeSum += kernel.size();
		if (kernel.size() > largestHardKernel) largestHardKernel = kernel.size();		
		
		// Start branching
		time = System.nanoTime();
		stage = 5;
		branch(1, vcrResult);
		branchTime += System.nanoTime() - time;
		stage = 6;

		if (curBranches > maxBranches) maxBranches = curBranches;
	}
	
	// Returns whether recursive deletion eliminated center vertex and kernel should be discarded
	private static boolean deleteKernelVertices (HashSet<Integer> deleteSet) {
		while (!deleteSet.isEmpty()) {
			int next = 0;								// dummy initialization
			for (int v : deleteSet) {next = v; break;}	// grab any vertex in delete set
			if (next == 0) return true;					// must delete center vertex
			deleteSet.remove(next);
			for (int n : kernel.getNeighbors(next)) {
				if (kernel.degree(n) < t || kernel.degree(n) < 2)
					deleteSet.add(n);
				if ((usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL) && n == 0)	// n is center vertes in kernel; 
					for (int n2 : kernel.getNeighbors(next))
						if (n2 != 0 && !needNoNeighbors(n, n2))
						{
							if (isInConflict(n, n2, false, 1)) // without 'next' 'n2' gets in conflict with center (only applies to HEREDITARY and VB_MODEL
								deleteSet.add(n2);
						}
			}
			kernel.deleteVertex(next);
		}
		return false;
	}
	
	private static void branch (int branchDepth, int lastVCRresult) {
		
		if (timeout != -1 && System.currentTimeMillis() - start > timeout * 1000) timeRunOut = true;
		if (timeRunOut) return;
		
		if (branchDepth > depth) depth++;
		branches++;
		curBranches++;
		
		// Abort branch if remaining graph is too small
		if (kernel.size() <= lowerBound) return;
		
		// Report solution if all vertices are marked
		if (marked.size() == kernel.size()) {
			if (usedModel == Model.HEREDITARY || marked.size() > t)
			{ 	// less than t vertices in the solution are only acceptable in the hereditary model ...
				solution = kernel.getClone();
				lowerBound = kernel.size();
			}
			return;
		}
		
		int minAgreeVertex = -1;
		// Determine vertex of minimum agreement
		for (int vertex : kernel.getVertices())
			if (!marked.contains(vertex)) {
				if (minAgreeVertex == -1 || agreements[vertex] < agreements[minAgreeVertex])
					minAgreeVertex = vertex;
		}
		
		// If vertex of minimum agreement equals size of kernel - 1, we found a solution (all vertices are agreeable)
		if (agreements[minAgreeVertex] == kernel.size() - 1 && (usedModel == Model.HEREDITARY || agreements[minAgreeVertex] > t-1)) {
			solution = kernel.getClone();
			lowerBound = kernel.size();
			return;
		}
		
		// Apply reduction rule R5 (Vertex Cover Rule) to see if branch can yield a better solution
		// (The VCR is relatively expensive, so we only want to check for it if absolutely necessary.)
		int currentVCRresult = lastVCRresult;
		// If no new conflicts were created since the last application of the VCR,
		// the last result of the VCR is still a valid upper bound for the number of required vertex deletions.
		// (Note that we only count conflicts as new which are not covered by the last VC of the conflict graph.)
		if (newConflicts == 0) {
			if (kernel.size() - lastVCRresult <= lowerBound) {
				// But some conflicts might have been resolved, so if lastVCRresult signals to abort the branch,
				// we need to run the VCR again to check whether the required number of vertex deletions has decreased
				currentVCRresult = vertexCoverRule(kernel.size() - lowerBound - 1);
				if (kernel.size() - currentVCRresult <= lowerBound) return;
			}	
		} else {
		// If X new conflicts were created, a new application of the VCR can at most report	lastVCRresult + X required
		// vertex deletions. We only need to check the VCR again, if lastVCRresult + X deletions shrink the kernel too much.
			if (kernel.size() - lastVCRresult - newConflicts <= lowerBound) {
				currentVCRresult = vertexCoverRule(kernel.size() - lowerBound - 1);
				newConflicts = 0;	// reset new conflicts counter
				if (kernel.size() - currentVCRresult <= lowerBound) return;
			}
		}
		
		// Use minimum agreement vertex as branch vertex
		int branchVertex = minAgreeVertex;
		
		// First branch over deleting the branch vertex
		int rollbackPoint = -1;
		if (!rollbackStack.isEmpty()) rollbackPoint = rollbackStack.peek();
		int newConflictsRestore = newConflicts;
		HashSet<Integer> lastVCRestore = lastVC;
		long time = System.nanoTime();
		boolean branch = delete(branchVertex);
		opsTime += System.nanoTime() - time;
		if (branch)
			branch(branchDepth+1, currentVCRresult);	// only branch if deleting created no unresolvable conflicts
		time = System.nanoTime();
		rollback(rollbackPoint);						// rollback changes of 1st child branch
		opsTime += System.nanoTime() - time;
		newConflicts = newConflictsRestore;				// also reset new conflicts counter to value before branch
		lastVC = lastVCRestore;							// as well as the last known vertex cover of the conflict graph
		// Then branch over marking the branch vertex
		time = System.nanoTime();
		branch = mark(branchVertex);
		opsTime += System.nanoTime() - time;
		if (branch) 
			branch(branchDepth+1, currentVCRresult);	// only branch if marking created no unresolvable conflicts
		time = System.nanoTime();
		rollback(rollbackPoint);						// rollback changes of 2nd child branch
		opsTime += System.nanoTime() - time;
		return;											// Done with this branch.
	}
	
	/**
	 * 
	 * @param vertex
	 * @param conflicts
	 * @return	Whether exploring this branch still makes sense.
	 * @throws AlgorithmException
	 * @throws GraphException
	 */
	private static boolean mark (int vertex) {
		
		ops++;
		marked.add(vertex);
		rollbackStack.push(vertex);
		// Apply reduction rule R3: Conflict resolution for newly marked vertex
		HashSet<Integer> recursiveDeleteSet = new HashSet<Integer>();
		for (int v : conflicts.get(vertex))
			if (marked.contains(v))
				return false;
			else if (kernel.contains(v))
				recursiveDeleteSet.add(v);				 		// Delete vertices in conflict with newly marked vertex
		if (kernel.size() - recursiveDeleteSet.size() <= lowerBound) return false;
		for (int v : recursiveDeleteSet)						// Perform recursive deletion
			if (kernel.contains(v))								// Recursive deletion might have already deleted v
				if (!delete(v)) return false;					// Abort branch if recursive deletion triggered R1
		if (kernel.size() <= lowerBound) return false;			// Check whether recursive deletion made kernel too small
		
		if (usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL)
		{
			// Apply reduction rule R6 (No Choice Rule)
			HashSet<Integer> mustMark = new HashSet<Integer>();
			for (int v : marked)
				if (v != vertex && !needNoNeighbors(vertex, v)) {
//					int minID = v<vertex?v:vertex;
//					int maxID = v>vertex?v:vertex;
					if (isInConflict(v, vertex, true, 1) && !isInConflict(v, vertex, true, 0)) //  noCommonNeighbors[maxID][minID] == t
						for (int must : kernel.getCommonNeighbors(vertex, v))
							if (!marked.contains(must))
								mustMark.add(must);
			}		
			for (int mustA : mustMark)								// First check must mark set for consistency
				for (int mustB : mustMark)
					if (mustA < mustB && !needNoNeighbors(mustA, mustB))
						if (isInConflict(mustB,mustA,true,0)) // noCommonNeighbors[mustB][mustA] < t
							return false;							// Cannot mark inconsistent set, abort branch
			for (int must : mustMark)								// Then check for consistency with existing marks
				for (int is : marked) {
//					int minID = must<is?must:is;
//					int maxID = must>is?must:is;
					if(isInConflict(must,is,true,0)) // !needNoNeighbors(must, is) && noCommonNeighbors[maxID][minID] < t
						return false;								// New marks lead to marked conflict, abort
				}
			for (int must : mustMark) {
				if (!kernel.contains(must))							// Other recursive markings might have forced its deletion
					return false;									// Must mark and can't, abort branch
				if (!marked.contains(must))							// Other recursive markings might have already marked it
					if (!mark(must))	return false;				// Abort branch if any recursive marking triggered R1
			}
		}
		completedOps++;
		return true;											// Continue exploring branch if R1 wasn't triggered
	}
	
	private static void unmark(int vertex) {
		marked.remove(vertex);
	}
	
	/**
	 * 
	 * @param vertex
	 * @return	Whether exploring this branch still makes sense.
	 * @throws AlgorithmException
	 * @throws GraphException
	 */
	private static boolean delete (int vertex) {
		
//		System.out.println(vertex);
		boolean carryOn = true;
		HashSet<Integer> recursiveDeleteSet = new HashSet<Integer>();
		HashSet<Integer> mustMark = new HashSet<Integer>();
		ops++;
		
		// Check whether deletion makes adjacent marked vertices unsupportable
		HashSet<Integer> adjacentMarked = new HashSet<Integer>();
		if (usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL)
		{
			for (int v : kernel.getNeighbors(vertex))
				if (marked.contains(v)) adjacentMarked.add(v);

			for (int v : adjacentMarked) 
				for (int w : adjacentMarked) 
					if (v < w && isInConflict(v, w, true, 1)) // !needNoNeighbors(v, w) && noCommonNeighbors[w][v] <= t
//						if (noCommonNeighbors[w][v] <= t)		// marked vertices will become unsupportable after deleting vertex
						return false;
		}
		
		// Perform reduction rule R2: Mark neighbors with too low degree for recursive deletion
		for (int i : kernel.getNeighbors(vertex))
			if (kernel.degree(i) <= t || kernel.degree(i) == 2)	{
				// will be t-1 (respectively 1) after deletion of original vertex
				if (marked.contains(i)) return false;	// cannot delete original vertex, abort branch
				recursiveDeleteSet.add(i);
			}
		
		if (kernel.size() - recursiveDeleteSet.size() <= lowerBound) return false;
		
		if (usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL)
		{// now adjust agreements-array for model HEREDITARY or VB_MODEL
			// Lower agreements of adjacent vertices by 1 (not for VB, they will be handled below)
			for (int v : kernel.getNeighbors(vertex)) {
				if (needNoNeighbors(v, vertex)) agreements[v]--;
				// Add vertex v to deletion set if it triggers reduction rule R4 (Purge vertices of too low agreement)
				if (agreements[v] < lowerBound) {
					if (marked.contains(v)) carryOn = false;	// Must purge marked vertex, abort branch
					else recursiveDeleteSet.add(v);
				}	
			}
			
			if (kernel.size() - recursiveDeleteSet.size() <= lowerBound)
				carryOn = false;								// Recursive delete set has become too large
					
			if (!carryOn) {									// If R4 triggered on marked vertex
				for (int v : kernel.getNeighbors(vertex))	// rollback changes so far and abort
					agreements[v]++;
				return false;								
			}
	
			
			// Lower agreements of non-adjacent (or all for VB) agreeable vertices by 1
			for (int v : kernel.getVertices())
				if (!needNoNeighbors(v, vertex) && v != vertex) {
//					int minID = v<vertex?v:vertex;
//					int maxID = v>vertex?v:vertex;
					if (!isInConflict(v, vertex, true, 0)) { //noCommonNeighbors[maxID][minID] >= t
						agreements[v]--;
						// Add vertex v to deletion set if it triggers reduction rule R4 (Purge vertices of too low agreement)
						if (agreements[v] < lowerBound) {				// Agreement became too low
							if (marked.contains(v))	carryOn = false;	// Must purge marked vertex, abort branch
							else recursiveDeleteSet.add(v);
						}
				}
			}		
			
			if (kernel.size() - recursiveDeleteSet.size() <= lowerBound)
				carryOn = false;								// Recursive delete set has become too large
					
			if (!carryOn) {									// If R4 triggered on marked vertex
				// rollback changes so far and abort
				if (usedModel == Model.HEREDITARY)							// we did not make these changes when using VB 
					for (int v : kernel.getNeighbors(vertex))
						agreements[v]++;						
				for (int v : kernel.getVertices())
					if (v != vertex && !needNoNeighbors(v, vertex)) {
//						int minID = v<vertex?v:vertex;
//						int maxID = v>vertex?v:vertex;
						if (!isInConflict(v, vertex, true, 0)) //noCommonNeighbors[maxID][minID] >= t
							agreements[v]++;
					}	
				return false;
			}
		
			// From now on we cannot abort mid-operation as noCommonNeighbors[][] might be left inconsistent
			for (int v : kernel.getNeighbors(vertex)) {
				for (int w : kernel.getNeighbors(vertex)) {				
					if (w <= v) continue;
					if (needNoNeighbors(v, w)) continue;
					noCommonNeighbors[w][v]--;
					if (isInConflict(v, w, true, 0) && !isInConflict(v, w, true, -1)) {							// Lost agreeability just now
//					if (noCommonNeighbors[w][v] == t-1) {							// Lost agreeability just now
						carryOn = notAgreeable(carryOn, recursiveDeleteSet, v, w);
					}
					
					// Collect vertices for reduction rule R6 (No Choice Rule)
					if (carryOn)		// Do nothing if branch is already marked for abortion
						if (adjacentMarked.contains(v) && adjacentMarked.contains(w) && isInConflict(v, w, true, 1) && !isInConflict(v, w, true, 0)) { // !needNoNeighbors(v, w) && noCommonNeighbors[w][v] == t
							for (int must : kernel.getCommonNeighbors(v, w)) 
								if (!marked.contains(must) && must != vertex)
									mustMark.add(must);
						}
				}	
			}
		}
		
		if (usedModel == Model.BICONNECTED)
		{
			for (int v : kernel.getVertices()) 
			{
				if (v != vertex && agree[Math.max(v, vertex)][Math.min(v, vertex)]) // v agrees with vertex but vertex will be removed...
					agreements[v]--;
			}
		}
		
		// We must perform deletion, even if we immediately undo it, to preserve consistency with noCommonNeighbors[][]
		kernel.deleteVertex(vertex);
		rollbackStack.push(vertex);
		
		if (usedModel == Model.BICONNECTED) 
		{
			for (int v : kernel.getNeighbors(vertex)) 
			{
				for (int w : kernel.getNeighbors(vertex)) 
				{
					if (w <= v) continue;
					noCommonNeighbors[w][v]--;
				}
			}
			for (int v : kernel.getVertices()) 
			{
				for (int w : kernel.getVertices()) 
				{
					if (w > v && agree[w][v]) // if w and v did not agree before deleting vertex, then they do not agree now...
					{
						connectednessLB[w][v]--; // possibly one less path between w and v
						connectednessReset[w][v] = false;
						if (isInConflict(w, v, true, 0))
						{
							disagreeReason[w][v] = vertex;
							carryOn = notAgreeable(carryOn, recursiveDeleteSet, v, w);
						}
					}
				}
			}
		}
		
		if (!carryOn) return false;
		// Check whether recursive delete set is too large
		if (kernel.size() - recursiveDeleteSet.size() <= lowerBound) return false;

		if (usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL)
		{
			// Check whether recursive marking and recursive deletion sets overlap -> R1 (Marked Conflict)
			for (int must : mustMark)
				for (int cant : recursiveDeleteSet)
					if (must == cant) return false;
			// Perform recursive deletion for reduction rules R2, R3 and R4
			for (int v : recursiveDeleteSet)						// Perform recursive deletion
				if (kernel.contains(v))								// Recursive deletion might have already deleted v
					if (!delete(v)) return false;					// Abort branch if recursive deletion triggered R1
			// Check whether recursive deletion deleted one of the vertices we must mark
			for (int must : mustMark)
				if (!kernel.contains(must))
					return false;									// Must mark deleted vertex, abort branch
			// Perform recursive marking according to reduction rule R6 (No Choice Rule)
			for (int mustA : mustMark)								// First check must mark set for consistency
				for (int mustB : mustMark)
					if (mustA < mustB)
						if (isInConflict(mustA, mustB, true, 0))
//						if (!needNoNeighbors(mustA, mustB) && noCommonNeighbors[mustB][mustA] < t)
							return false;							// Cannot mark inconsistent set, abort branch
			for (int must : mustMark)								// Then check for consistency with existing marks
			{
				for (int is : marked) if (must != is) {
//					int minID = must<is?must:is;
//					int maxID = must>is?must:is;
					if (isInConflict(must, is, true, 0))
//					if(!needNoNeighbors(must, is) && noCommonNeighbors[maxID][minID] < t)
						return false;								// New marks lead to marked conflict, abort
				}
			}
			for (int must : mustMark) {
				if (!kernel.contains(must))				// Other recursive markings might have forced its deletion
					return false;						// Must mark and can't, abort branch	
				if (!marked.contains(must))				// Other recursive markings might have already marked it
					if (!mark(must))					// Abort branch if any recursive marking triggered R1
						return false;					// Must mark and can't, abort branch
			}
		}
		completedOps++;
		return true;								// Continue exploring branch if R1 wasn't triggered
	}

	/**
	 * adjusts agreement values if v and w (w > v !!) are now in conflict
	 * @param carryOn
	 * @param recursiveDeleteSet
	 * @param v
	 * @param w
	 * @return
	 */
	private static boolean notAgreeable(boolean carryOn, HashSet<Integer> recursiveDeleteSet, int v, int w) 
	{
		agree[w][v] = false;
		agreements[v]--;
		agreements[w]--;
		if (!carryOn) return carryOn;					// Everything below this point is not needed for rollback
		conflicts.get(v).add(w);						// Deleting vertex generates conflict between v and w
		conflicts.get(w).add(v);
		if (marked.contains(v)) recursiveDeleteSet.add(w);		// Delete vertices in conflict with marked ones
		else if (marked.contains(w)) recursiveDeleteSet.add(v);
		if (lastVC != null)	
			if (!lastVC.contains(v) && !lastVC.contains(w))	// if our last VC of the conflict graph does not cover
				newConflicts++;								// the new conflict increase new conflicts counter by 1
		// Collect vertices for reduction rule R4 (Purge vertices of too low agreement)
		if (agreements[v] < lowerBound) {				// Agreement became too low
			if (marked.contains(v)) carryOn = false;	// Must purge marked vertex, abort branch
			else recursiveDeleteSet.add(v);
		}	
		if (agreements[w] < lowerBound) {				// Agreement became too low
			if (marked.contains(w)) carryOn = false;	// Must purge marked vertex, abort branch
			else recursiveDeleteSet.add(w);
		}
		return carryOn;
	}

	
	private static void undelete (int vertex) {
		
		kernel.undeleteVertex(vertex);
		if (usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL)
		{
			// Increase agreements of non-adjacent agreeable (or all for VB) vertices by 1
			for (int v : kernel.getVertices())
				if (v != vertex && !needNoNeighbors(v, vertex)) {
//					int minID = v<vertex?v:vertex;
//					int maxID = v>vertex?v:vertex;
					if (!isInConflict(v, vertex, true, 0)) //noCommonNeighbors[maxID][minID] >= t
						agreements[v]++;
			}		
					
			for (int v : kernel.getNeighbors(vertex)) {
				if (needNoNeighbors(v, vertex))
					agreements[v]++;							// Increase agreement of adjacent vertices by 1
				for (int w : kernel.getNeighbors(vertex)) {
					if (w <= v) continue;						// Avoids considering pair (v,w) twice
					if (needNoNeighbors(v, w)) continue;		// No need to consider adjacent vertices
					noCommonNeighbors[w][v]++;
					if (isInConflict(v, w, true, 1) && !isInConflict(v, w, true, 0)) {
//					if (noCommonNeighbors[w][v] == t) {						// Gained agreeability just now
						conflicts.get(v).remove(w);
						conflicts.get(w).remove(v);
						agreements[v]++;
						agreements[w]++;
					}
				}	
			}
		}
		if (usedModel == Model.BICONNECTED) 
		{
			for (int v : kernel.getNeighbors(vertex)) 
			{
				for (int w : kernel.getNeighbors(vertex)) 
				{				
					if (w <= v) continue;
					noCommonNeighbors[w][v]++;
				}
			}
			for (int v : kernel.getVertices()) 
			{
				if (v != vertex && agree[Math.max(v, vertex)][Math.min(v, vertex)]) // v agreed with vertex but vertex was removed...
					agreements[v]++;
				
				for (int w : kernel.getVertices()) 
				{
					if (w > v)
					{
						if (!connectednessReset[w][v])
						{
							connectednessLB[w][v]++;
							connectednessReset[w][v] = false;
						}
						if (disagreeReason[w][v] == vertex) // if w and v did not agree after deleting vertex (before they agreed))
						{
							conflicts.get(v).remove(w);
							conflicts.get(w).remove(v);
							agreements[v]++;
							agreements[w]++;
							agree[w][v] = true;
							disagreeReason[w][v] = -1;
						}
					}	
				}
			}
		}
	}

	private static int vertexCoverRule (int limit) {
	    if (useTriangleBound) 
	    	return trianglePackingRule(limit);
	    else 
	    	return matchingRule(limit);
	}

	private static int matchingRule (int limit) {
		vcrApplications++;
		long time = System.nanoTime();
		// 2-approximation for vertex cover of conflict graph
		// Minimal vertex cover of conflict graph is lower bound on vertices to delete 
		HashSet<Integer> marks = new HashSet<Integer>();
		int size = 0;
		for (int v : conflicts.keySet())	
			if (kernel.contains(v) && !marks.contains(v))	
				for (int w : conflicts.get(v))		
					if (kernel.contains(w) && w > v)					// every remaining conflict is an edge
						if (!marks.contains(v) && !marks.contains(w)) {	// if the edge is not covered	
							marks.add(v); marks.add(w); size++;			// cover edge with both end points
							if (size > limit) return size;				// if current size exceeds allowed limit
							break;										// v us marked, so the condition of the if-statement is never met again.
						}												// of deletions we do not need to carry on
		lastVC = marks;													// remember this vertex cover
		vcrTime += System.nanoTime() - time;
		return size;															// half the size of the 2-approximation
	}

   private static int trianglePackingRule (int limit) 
   {
	   vcrApplications++;
	   long time = System.nanoTime();
	   // Compute Lower bound on vertices to delete: Triangles need 2 deletions, edges need 1 deletion 
	   HashSet<Integer> marks = new HashSet<Integer>();
	   int size = 0;
	   //First checking for conflict triangles
	   for (int v : conflicts.keySet()) 
		   if (kernel.contains(v) && !marks.contains(v))	
			   for (int w : conflicts.get(v))
			   {
				   if(marks.contains(v))
				   {
					   break; // v us marked, so the condition of the next if-statement is never met again.
				   }
				   if (kernel.contains(w) && w > v && !marks.contains(w) && !marks.contains(v))	    // w is in conflict with and unmarked v
				   {
				       HashSet<Integer> triangleCand =  new HashSet<Integer>(conflicts.get(w)); 
				       triangleCand.retainAll(conflicts.get(v));
				       for (int x : triangleCand)  									//x is in conflict with w and v
					   if (kernel.contains(x) && x > w)
					       if (!marks.contains(x) && !marks.contains(w) && !marks.contains(v)) {// if the edge is not covered	
							   marks.add(v); marks.add(w); marks.add(x); size +=2;	// cover edge with both end points
							   if (size > limit) return size;						// if current size exceeds allowed limit
							   break; 												// v us marked, so the condition of the if-statement is never met again.
					       }									        			// of deletions we do not need to carry on
				   }
			   }
		//Now checking for conflict pairs
		for (int v : conflicts.keySet())	
			if (kernel.contains(v) && !marks.contains(v))	
				for (int w : conflicts.get(v))		
					if (kernel.contains(w) && w > v)	// every remaining conflict is an edge
						if (!marks.contains(v) && !marks.contains(w)) {		// if the edge is not covered		
							marks.add(v); marks.add(w); size++;				// cover edge with both end points
							if (size > limit) return size;					// if current size exceeds allowed limit
							break;											// v us marked, so the condition of the if-statement is never met again.
						}													// of deletions we do not need to carry on
		lastVC = marks;														// remember this vertex cover
		vcrTime += System.nanoTime() - time;
		return size;														// half the size of the 2-approximation
	}

	private static void rollback (int rollbackPoint) {
		while (!rollbackStack.empty() && rollbackStack.peek() != rollbackPoint) {
			int vertex = rollbackStack.pop();
			if (marked.contains(vertex)) unmark(vertex);	// If vertex is marked operation to undo must be marking
			else undelete(vertex);							// otherwise it must be a deletion
		}	
	}
	
	// Helper method to distinguish between the models, always returns false
	// for Veremyev/Boginski and returns false for our model iff vertices are not adjacent. 
	private static boolean needNoNeighbors (int vertexA, int vertexB) {
		return (usedModel == Model.HEREDITARY && kernel.adjacent(vertexA, vertexB)); 
	}
	
	/**
	 * checks if two vertices in the graph kernel are in conflict wrt. the used model. 
	 * 
	 * For our model (Model.HEREDITARY) the vertices vertexA and vertexB are in conflict if they are not adjacent and have less than t vertices.
	 * 
	 * For the VB model (Model.VB_MODEL) the vertices vertexA and vertexB are in conflict if there are less than t vertex-disjoint paths of length two between them, that is, 
	 * vertexA and vertexB have less than t-1 common neighbors (if vertexA and vertexB are adjacent) or less than t common neighbors (if vertexA and vertexB are not adjacent).
	 * 
	 * For the hereditary model (usedModel == Model.BICONNECTED) the vertices vertexA and vertexB are in conflict if they are at distance more than 2 (for the other two 
	 * models this is checked implicitly by the common neighbors) or if there are less than t vertex-disjoint paths (of any length) between them.

	 * @param vertexA
	 * @param vertexB
	 * @param useCommonNeighborArray Can the array noCommonNeighbors be used? Set to false if the array is not yet build or the graph changed
	 * @param deleteCommonNeighbors how many common neighbors of vertexA and vertexB will be removed
	 * 			
	 * @return
	 */
	private static boolean isInConflict (int vertexA, int vertexB, boolean useCommonNeighborArray, int deleteCommonNeighbors) {
		if (vertexA == vertexB)
			return false;
		conflictChecks++;
		
		int maxId = Math.max(vertexA, vertexB);
		int minId = Math.min(vertexA, vertexB);
		
		int noOfNeighbors = 0;
		if (useCommonNeighborArray)
			noOfNeighbors = noCommonNeighbors[maxId][minId];
		else
			noOfNeighbors = kernel.countCommonNeighbors(vertexA, vertexB);
		noOfNeighbors = noOfNeighbors - deleteCommonNeighbors;
		
		if (usedModel == Model.HEREDITARY)
		{
			return !kernel.adjacent(vertexA, vertexB) && noOfNeighbors < t;
		} 
		else if (usedModel == Model.VB_MODEL)
		{
			if (kernel.adjacent(vertexA, vertexB))
				return noOfNeighbors < t-1;
			else
				return noOfNeighbors < t;
		} 
		else if (usedModel == Model.BICONNECTED)
		{
			if (!kernel.adjacent(vertexA, vertexB) && noOfNeighbors <= 0)
				return true; // vertices at distance more than two => they are in conflict
			
			if (noOfNeighbors >= t)
				return false; // vertices have at least t common neighbors and thus at least t internally vertex disjoint paths => not in conflict
			
			boolean result;
			long time = System.nanoTime();
			
			if (!useCommonNeighborArray)
			{
				result = !kernel.enoughInternallyVertexDisjointPaths(vertexA, vertexB, t);
				maxFlowCallsPerStage[stage]++;
				maxFlowCalls++;
			}
			else if (connectednessLB[maxId][minId] < t)
			{
				connectednessLB[maxId][minId] = kernel.noOfInternallyVertexDisjointPaths(vertexA, vertexB);
				connectednessReset[maxId][minId] = true;
				result = connectednessLB[maxId][minId] < t;
				maxFlowCallsPerStage[stage]++;
				maxFlowCalls++;
			}
			else
			{
				return false;
			}
			timeSpentOnMaxFlow += System.nanoTime() - time;
			
			if (result)
				negAnswer++;
			else
				posAnswer++;
			
			return result;
		}
		return false;
	}
	
	/**
	 * function checking whether the array agreements is updated consistently. 
	 * @return
	 */
	private static boolean checkAgreements ()
	{
//		if (usedModel == Model.BICONNECTED)
		{
			int[] tmpAgreements = new int[agreements.length];
			for (int v : kernel.getVertices())
				for (int w : kernel.getVertices()) if (v < w)
					if  (!isInConflict(w, v, false, 0)) { 
						tmpAgreements[v]++;
						tmpAgreements[w]++;
					}
			for (int v : kernel.getVertices())
			{
				if (tmpAgreements[v] != agreements[v])
				{
					System.out.println("Agreement incorrect?! Check vertex " +v);
					System.out.println("tmpAgreements: "+ arrayToString(tmpAgreements));
					System.out.println("agreements:    "+ arrayToString(agreements));
					System.out.println("vertices: "+ kernel.getVertices());
//					System.exit(0);
					return false;
				}
			}
		}
		return true;
	}
	
	private static int tmpCounter = 0;
	private static String arrayToString (int[] a)
	{
		tmpCounter++;
        if (a == null)
            return "null";
        int iMax = a.length - 1;
        if (iMax == -1)
            return "[]";

        StringBuilder b = new StringBuilder();
        b.append('[');
        for (int i = 0; ; i++) {
            b.append(i +": " + a[i]);
            if (i == iMax)
                return b.append("] " + tmpCounter).toString();
            b.append(", ");
        }
	}
}
