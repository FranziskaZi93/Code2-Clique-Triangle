package robustTwoClub.correctness;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import robustTwoClub.algorithms.RobustTwoClubAlgorithm;
import robustTwoClub.algorithms.RobustTwoClubAlgorithm.Model;
import robustTwoClub.graph.RtcGraph;



/** This class implements an inefficient solver for the 2,t-CLUB problem, which is simple and easily
 *  checkable for correctness. Its purpose is to detect errors in the actual implementation
 *  of the 2,t-CLUB solver by performing an exhaustive search for solutions larger than
 *  the reported optimum solution on small test instances.  
 *
 */
public class CorrectnessChecker {
	
	private static RtcGraph graph;
		
	private static boolean check2tClub (RtcGraph club, int t, String algoName, Model usedModel, boolean printMessage) {
		
		int hubvertex = -1;
		int noHubs = 0;
		for (int v : club.getVertices())
		{
			for (int w : club.getVertices())
				if (v != w && (!club.adjacent(v, w) || usedModel != Model.HEREDITARY)) {
					int countCommonNeighbors = club.countCommonNeighbors(v, w);
					if (club.adjacent(v, w))
						countCommonNeighbors++;
					
					if (countCommonNeighbors == 0)
					{
						if (printMessage)
						{
							System.out.println("Checker says: Solution of algorithm is INCORRECT!");
							System.out.println("Algorithm "+algoName+" claimed this to be a 2,"+t+"-club: " + club.getVertexNames().toString());
							System.out.println("But this is not true. Vertices "+club.getVertexName(v)+" and "+ club.getVertexName(w) + " are at distance more than 2!");
							return false;
						}
						return false;
					}
					
					if ( (countCommonNeighbors < t && (usedModel == Model.VB_MODEL || usedModel == Model.HEREDITARY))
							|| (usedModel == Model.BICONNECTED && !club.enoughInternallyVertexDisjointPaths(v, w, t)) )
					{
						if (printMessage)
						{
							System.out.println("Checker says: Solution of algorithm is INCORRECT!");
							System.out.println("Algorithm "+algoName+" claimed this to be a 2,"+t+"-club: " + club.getVertexNames().toString());
							System.out.println("But this is not true. Vertices "+club.getVertexName(v)+" and "+ club.getVertexName(w) + " violate the 2,"+t+"-club property!");
							return false;
						}
						return false;
					}
				}
			if (usedModel == Model.BICONNECTED && t == 2 && club.degree(v) + 1 == club.size()) //extra check
			{
				hubvertex = v;
				noHubs++;
			}
		}
		if (usedModel == Model.BICONNECTED && noHubs == 1 && t == 2)
		{
			RtcGraph copy = club.getClone();
			copy.deleteVertex(hubvertex);
			if (!copy.isConnected())
			{
				if (printMessage)
				{
					System.out.println("Checker says: Solution of algorithm is INCORRECT!");
					System.out.println("Algorithm "+algoName+" claimed this to be a biconnected 2,"+t+"-club: "+ club.getVertexNames().toString()+"  but removing " + hubvertex + " disconnects the graph.");
				}
				return false;
			}		
		}
		if (usedModel != Model.HEREDITARY && club.size() <= t)
		{
			if (printMessage)
			{
				System.out.println("Checker says: Solution of algorithm is INCORRECT!");
				System.out.println("Algorithm "+algoName+" claimed this to be a biconnected 2,"+t+"-club: "+ club.getVertexNames().toString()+"  but it has size " + club.size());
			}
			return false;
		}
		return true;				
	}
	
	/** Checks whether the claim "Graph g has no 2,t-club of size > opt" is true by performing
	 *  an computationally inefficient (but easily verifiable) brute force search for a better
	 *  solution. Only feasible for small instances (of about 20-30 vertices, depending on density).
	 *  Outputs following information to the console, if claim is not true:
	 *  	1) Size of one 2,t-club in g larger than opt. (not necessarily the largest solution)
	 *  	2) Vertices contained in the found better solution.
	 *  Checker only works on graphs with vertex IDs 0 to N-1, where N is the number of vertices.
	 *  
	 * @param g			Graph g of the 2,t-CLUB instance.
	 * @param t			Parameter t of the 2,t-CLUB instance.
	 * @param opt		Claimed size of best solution for 2,t-CLUB on graph g.
	 * @return			Whether the claim is true.
	 * @throws 			GraphException
	 * @throws 			AlgorithmException
	 */
	private static boolean testOptimality (RtcGraph g, int t, int opt, Model usedModel) {
		graph = g.getClone();			// Make a local clone in order to leave original graph unchanged
		HashSet<Integer> marks = new HashSet<Integer>();
		
		// Assert that graph has vertex IDs 0 to N-1
		for (int i = 0; i < g.size(); i++)
			if (!g.contains(i)) {
				System.out.println("STOP! Checker only works on graphs with integer vertex IDs from 0 to N-1.");
				return false;
			}	
		
		for (int size = opt + 1; size <= g.size(); size++)
			if (recursiveMarking(size, marks, -1, t, usedModel)) break;	
			// Abort exhaustive search if we find better solution
		if (marks.isEmpty()) {
			System.out.println("Checker says: Solution of algorithm is CORRECT.");
			return true;
		} else {
			System.out.println("Checker says: Solution of algorithm is INCORRECT!");
			System.out.println("Algorithm claimed there is no solution of size > "+opt+", ");
			System.out.print("but checker found following size "+marks.size()+" solution: ");
			System.out.println(marks.toString());
			return false;
		}
	}
	
	/** Recursively marks vertices until desired number of marks has been placed,
	 *  then calls testMarks() to determine whether marked set is a 2,t-CLUB solution.
	 *  
	 * @param M			Number of marks left.
	 * @param marks		Vertices marked so far.
	 * @param lastMark	Index of last mark added. Checker algorithm only adds marks in ascending index order,
	 * 					in order to avoid generating more than one permutation leading to the same subset.
	 * 					Call with lastMark = -1 as root for the recursion tree.
	 * @return			Whether some size M subset of marks together with old marks constituted a 2,t-club.
	 * @throws 			GraphException 
	 */
	private static boolean recursiveMarking (int M, HashSet<Integer> marks, int lastMark, int t, Model usedModel) {
		
		if (M == 0)	return testMarks(marks, t, usedModel);	// If all marks are set, check for 2,t-clubness
		
		// Recursively generate all subsets of M vertices, such that ...
		for (int newMark = lastMark+1; newMark < graph.size(); newMark++) {	// ... no vertex is selected twice, ... 		 
			boolean agreeable = true;
//			for (int oldMark : marks)										// ... and for every pair ...	
//				if (!graph.adjacent(oldMark, newMark) || usedModel == Model.VB_MODEL) {
//					int countCommonNeighbors = graph.countCommonNeighbors(oldMark, newMark);
//					if (graph.adjacent(oldMark, newMark) && usedModel == Model.VB_MODEL)
//						countCommonNeighbors++;
//					if (usedModel == Model.BICONNECTED)
//						countCommonNeighbors = countCommonNeighbors*t; // in this model one common neighbor is sufficient
//					if (countCommonNeighbors < t)
//						agreeable = false;
//				}									// ... there are at least t common neighbors in G
			if (agreeable) {
				marks.add(newMark);
				// Abort exhaustive search as soon as we find some solution
				if (recursiveMarking(M-1, marks, newMark, t, usedModel)) return true;						
				marks.remove(newMark);
			}
		}	
		return false;									// Exhaustive search for subset 'marks' yielded no solution
	}
	
	/** Checks whether a given set of marked vertices is a 2,t-club.
	 * 
	 * @param marks		Set of marked vertices.
	 * @return			Whether set of marked vertices constitutes a 2,t-club.
	 * @throws 			GraphException
	 */
	private static boolean testMarks (HashSet<Integer> marks, int t, Model usedModel) {
		RtcGraph subset = graph.getClone();
		for (int vertex : graph.getVertices())
			if (!marks.contains(vertex))
				subset.deleteVertex(vertex);
		return check2tClub(subset, t, "", usedModel, false);
	}
	
	/** Runs multiple tests of the algorithm on randomly generated Erdös-Renyi graphs.
	 *  For each test a random graph with 'minN' - 'maxN' vertices and random edge
	 *  probability p with minP <= p <= maxP is generated and the 2,t-CLUB algorithm
	 *  is invoked to solve 2,t-CLUB for random t with minT <= t <= maxT.
	 *  The claimed optimum will then be checked by the optimum checker.
	 * 
	 * @param tests		Number of tests to run.
	 * @param minN		Minimum size (# of vertices) of generated graphs.
	 * @param maxN		Maximum size (# of vertices) of generated graphs.
	 * @param minP		Lower bound for edge probability of generated graphs.
	 * @param maxP		Upper bound for edge probability of generated graphs.
	 * @param minT		Lower bound for parameter t of problem.
	 * @param maxT		Upper bound for parameter t of problem.
	 * @param output	Whether graphs for which the algorithm failed should be written to a file.
	 * 					(Filename: "ErdösRenyiTest" + test number + ".edges")
	 * @param usedModel		Whether the Veremyev/Boginski model is used or not.
	 * @return			Whether all test were passed.
	 * @throws IOException 
	 * @throws GraphException 
	 * @throws AlgorithmException 
	 */
	public static boolean runErdosRenyiTests (int tests, int minN, int maxN, double minP, double maxP,
			int minT, int maxT, boolean output, Model usedModel) {
		boolean passedAll = true;
		
		PrintWriter out = null;
		try {
			FileWriter fw = new FileWriter("tests/errors.txt", true);
			BufferedWriter bw = new BufferedWriter(fw);
		    out = new PrintWriter(bw);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
		
		for (int test = 1; test <= tests; test++) {
//			int n = (int) Math.floor(Math.random() * (maxN - minN + 1) + minN);
			int n = (int) Math.floor(minN + (test-1)*(maxN - minN + 1)/tests);
			double p = Math.random() * (maxP - minP) + minP;
			int t = (int) Math.floor(Math.random() * (maxT - minT + 1) + minT);
			RtcGraph original = RtcGraph.ErdosRenyi(n, p);
			original.writeToFile("tests", "tmp", new String[1]);
			original = new RtcGraph("tests/tmp.dimacs",0); // writing to file and then reading from it removes degree-zero vertices from the graph...
			if (original.size() == 0) 
			{
				continue; //corner case our algorithm cannot deal with ;)
			}
//			original = new RtcGraph("tests/BICONNECTEDErdösRenyiTest42.dimacs",0); 
//			original = new RtcGraph("graphs/clustering/karate.graph.dimacs",0);
			RtcGraph algoCopy = original.getClone();
			RtcGraph algoNoKernelsCopy = original.getClone();
			RtcGraph checkerCopy = original.getClone();
			System.out.println("TEST  "+test+":\t Graph with "+original.size()+" vertices and "+original.getEdgeCount() +
					" edges (p = "+p+").   ");
			System.out.println("Task:\t\t Find size of largest 2,"+t+"-club ("+usedModel+").\n");
			RtcGraph solution = RobustTwoClubAlgorithm.run(algoCopy, "", t, t, 0, 0, 0, false, null, false, 0, -1, true, usedModel);
			RtcGraph solutionNK = RobustTwoClubAlgorithm.run(algoNoKernelsCopy, "", t, t, 0, 0, 0, false, null, true, 0, -1, true, usedModel);
			boolean passedTest = true;
			int solSize = 0, solSizeNK = 0;
			if ((solution == null && solutionNK != null) || (solution != null && solutionNK == null) ||
					(solution != null && solutionNK != null && (solution.size() != solutionNK.size()))) {
				if (solution == null) {
					System.out.println("Kernelization variant reported no solution, "+
								"so solution must be clique of size < "+(t+1)+".\n");
					passedTest = passedTest && testOptimality(checkerCopy, t, t, usedModel);
				} else {
					System.out.println("Kernelization variant reported optimal solution of size "+solution.size()+".");
					System.out.println("Vertices in solution: "+solution.getVertexNames().toString()+"\n");
					passedTest = passedTest && testOptimality(checkerCopy, t, solution.size(), usedModel);
					passedTest = passedTest && check2tClub(solution, t, "K", usedModel, true);
					solSize = solution.size();
				}
				if (solutionNK == null) {
					System.out.println("\nSingle pass variant reported no solution, "+
								"so solution must be clique of size < "+(t+1)+".\n");
					passedTest = passedTest && testOptimality(checkerCopy, t, t, usedModel);
				} else {
					System.out.println("\nSingle pass variant reported optimal solution of size "+solutionNK.size()+".");
					System.out.println("Vertices in solution: "+solutionNK.getVertexNames().toString()+"\n");
					passedTest = passedTest && testOptimality(checkerCopy, t, solutionNK.size(), usedModel);
					passedTest = passedTest && check2tClub(solutionNK, t, "NK", usedModel, true);
					solSizeNK = solutionNK.size();
				}
			} else {
				if (solution != null) passedTest = passedTest && check2tClub(solution, t, "K", usedModel, true);
				if (solutionNK != null) passedTest = passedTest && check2tClub(solutionNK, t, "NK", usedModel, true);
				if (passedTest) {
					if (solution == null) {
						System.out.println("Algorithm reported no solution, so solution must be clique of size < "+(t+1)+".\n");
						passedTest = testOptimality(checkerCopy, t, t, usedModel);
					} else {
						System.out.println("Algorithm reported optimal solution of size "+solution.size()+".");
						System.out.println("Vertices in solution: "+solution.getVertexNames().toString()+"\n");
						passedTest = testOptimality(checkerCopy, t, solution.size(), usedModel);
					}
				}
			}
			passedAll = passedAll && passedTest;
			if (!passedTest && output) {
				original.writeToFile("tests", usedModel+"ErdösRenyiTest"+test, new String[1]);
				out.println("Output graph on which algorithm failed to 'ErdösRenyiTest"+test+".dimacs'. \n");
				System.out.println("Output graph on which algorithm failed to 'ErdösRenyiTest"+test+".dimacs'.");
				out.println("Solution (kernel)   of size "+solSize+" when loading: " + (0<solSize?solution.getVertices():""));
				out.println("Solution (nokernel) of size "+solSizeNK+" when loading: " + (0<solSizeNK?solutionNK.getVertices():""));
				System.exit(0);
			}
			System.out.println("\n");
		}
		out.close();
		if (passedAll)
			System.out.println("Algorithm PASSED all tests.");
		else System.out.println("Algorithm FAILED tests.");
		return passedAll;
	}
	
	public static boolean clusteringInstanceTest() 
	{
		boolean passedAll = true;
		String prefix = "graphs/clustering/";
		String postfix = ".graph.dimacs";
		int[] tValues = {1,2,3,4,5,7,9,10,15,20,50,100,1000};
		HashMap<String, HashMap<Model, int[]>> allSolutions = new HashMap<String, HashMap<Model, int[]>>();
		HashMap<Model, int[]> solutions;
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {32,12,5,4,4,4,4,4,4,4,4,4,4});
		solutions.put(Model.BICONNECTED, new int[] {32,30,11,0,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {32,12,5,0,0,0,0,0,0,0,0,0,0});
		allSolutions.put("add32", solutions);

		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {50,23,12,9,5,5,5,5,5,5,5,5,5});
		solutions.put(Model.BICONNECTED, new int[] {50,48,44,39,33,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {50,23,12,6,0,0,0,0,0,0,0,0,0});
		allSolutions.put("adjnoun", solutions);

		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {6,4,3,3,3,3,3,3,3,3,3,3,3});
		solutions.put(Model.BICONNECTED, new int[] {6,6,0,0,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {6,3,0,0,0,0,0,0,0,0,0,0,0});
		allSolutions.put("cs4", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {13,9,7,6,5,5,5,5,5,5,5,5,5});
		solutions.put(Model.BICONNECTED, new int[] {13,12,12,11,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {13,9,7,6,0,0,0,0,0,0,0,0,0});
		allSolutions.put("dolphins", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {16,14,13,13,12,11,9,9,9,9,9,9,9});
		solutions.put(Model.BICONNECTED, new int[] {16,16,15,15,15,13,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {16,14,13,12,12,10,0,0,0,0,0,0,0});
		allSolutions.put("football", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {51,33,24,24,24,24,24,24,24,24,24,24,24});
		solutions.put(Model.BICONNECTED, new int[] {51,45,40,31,28,24,24,24,24,24,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {51,33,24,24,24,24,24,24,24,24,0,0,0});
		allSolutions.put("hep-th", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {18,12,6,6,5,5,5,5,5,5,5,5,5});
		solutions.put(Model.BICONNECTED, new int[] {18,17,12,9,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {18,12,6,6,0,0,0,0,0,0,0,0,0});
		allSolutions.put("karate", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {35,22,21,20,20,20,20,20,20,20,20,20,20});
		solutions.put(Model.BICONNECTED, new int[] {35,25,21,20,20,20,20,20,20,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {35,22,21,20,20,20,20,20,20,0,0,0,0});
		allSolutions.put("netscience", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {28,20,15,13,11,6,6,6,6,6,6,6,6});
		solutions.put(Model.BICONNECTED, new int[] {28,28,28,27,26,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {28,20,15,12,10,0,0,0,0,0,0,0,0});
		allSolutions.put("polbooks", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {20,9,7,6,6,6,6,6,6,6,6,6,6});
		solutions.put(Model.BICONNECTED, new int[] {20,14,12,11,6,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {20,9,7,6,6,0,0,0,0,0,0,0,0});
		allSolutions.put("power", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {5,4,3,3,3,3,3,3,3,3,3,3,3});
		solutions.put(Model.BICONNECTED, new int[] {5,5,0,0,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {5,3,0,0,0,0,0,0,0,0,0,0,0});
		allSolutions.put("uk", solutions);
		
		solutions = new HashMap<Model, int[]>();
		solutions.put(Model.HEREDITARY, new int[] {9,6,5,3,3,3,3,3,3,3,3,3,3});
		solutions.put(Model.BICONNECTED, new int[] {9,9,9,0,0,0,0,0,0,0,0,0,0});
		solutions.put(Model.VB_MODEL, new int[] {9,6,0,0,0,0,0,0,0,0,0,0,0});
		allSolutions.put("whitaker3", solutions);
//		
//		solutions = new HashMap<Model, int[]>();
//		solutions.put(Model.HEREDITARY, new int[] {});
//		solutions.put(Model.BICONNECTED, new int[] {});
//		solutions.put(Model.VB_MODEL, new int[] {});
//		allSolutions.put("", solutions);

		for (HashMap.Entry<String, HashMap<Model, int[]>> entry : allSolutions.entrySet())
		{
			String file = prefix + entry.getKey() + postfix;
			System.out.println("Testin File " + file);
			for (HashMap.Entry<Model, int[]> entry2 : entry.getValue().entrySet())
			{
				for (int i = 0; i < tValues.length; i++)
				{
					RtcGraph graph = new RtcGraph(file,0);
					RtcGraph solution = RobustTwoClubAlgorithm.run(graph, "", tValues[i], tValues[i], 0, 0, 0, false, null, false, 0, -1, true, entry2.getKey());
					int solSize = 0;
					if (solution != null)
						solSize = solution.size();
					if (solSize != entry2.getValue()[i])
					{
						passedAll = false;
						System.out.println("ERROR: graph: " + file + " \t Model: " + entry2.getKey() + " \t t: " + tValues[i]);
					}
				}
			}
		}
		if (passedAll)
			System.out.println("Algorithm PASSED all tests.");
		else 
			System.out.println("Algorithm FAILED tests.");
		return passedAll;
	}
}
