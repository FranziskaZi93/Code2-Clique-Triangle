package robustTwoClub.analysis;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import robustTwoClub.algorithms.CliqueAlgorithm;
import robustTwoClub.algorithms.RobustTwoClubAlgorithm.Model;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.latex.TableMaker;

public class RtcGraphAnalysis {

	public static String[] header;
	public static HashSet<Integer> clique;
	
	private static boolean timedOut;	// Timeout flag for algorithms for hard problems
	private static long start;			// Start time variable for timing algorithms for hard problems
	
	private static final DecimalFormat format5 = new DecimalFormat("0.00000");
	private static final DecimalFormat format5e2 = new DecimalFormat("0.00000E00");
	
	private static RtcGraph compareSolution;
	
	/** Analyzes given input graph. Outputs results to console, if 'silent' is not set true.
	 *  Additionally generates header information as string array for output files.
	 *  Returns size of maximum clique of graph, if analyzeLevel = 2.
	 * @param graph		The graph to be analyzed.
	 * @param level		Analysis level. 1 = Only basic information computable in linear time.
	 * 									2 = Additional information requiring higher order polynomial time.
	 * 									3 = Even more information requiring to solve NP-hard problems.
	 * @param name		Name of the graph to be analyzed. Usually pass file name here.
	 * @param silent	Whether output to the console shall be suppressed.
	 * @param timeout	Timeout in seconds for algorithms of NP-hard problems used by graph analysis.
	 * @param subCall	Set to true, if method was called by analyzeSolution.
	 * @return			Whether largest found clique of graph is guaranteed to be maximum.
	 * 					Header and largest clique can be retrieved from public variables 'header' and 'clique'
	 * 					after call to this method.
	 */
	public static boolean analyzeGraph (RtcGraph graph, int level, String name, boolean silent, int timeout, boolean subCall) {
		boolean maxClique = false;
		ArrayList<String> headerInfo = new ArrayList<String>();
		if (!subCall) headerInfo.add("Statistics for the input graph \"" + name + "\":");
		else headerInfo.add("Statistics for this solution: ");
		// Level 1 statistics
		long pairs = (long)(graph.size() - 1) * graph.size() / 2;
		int edgeCount = graph.getEdgeCount();
		double density = (double)(edgeCount) / pairs;
		String d = format5.format(density);
		if (density < 0.0001) d = format5e2.format(density);
		headerInfo.add("vertices: " + graph.size() + "   edges: " + edgeCount + 
				"   edge density: " + d);
		int degreeInfo[] = getDegreeInfo(graph);
		headerInfo.add("min. degree: " + degreeInfo[0] + "   max. degree: " + degreeInfo[1]
				+ "   avg. degree: " + degreeInfo[2] + "   h-index: " + degreeInfo[3]);
		// Level 2 statistics
		if (level > 1) {
			if (!subCall) {		// information pointless for solutions, which are all diameter 2 graphs
				int nbInfo[] = get2NeighborhoodInfo(graph);
				headerInfo.add("2-neighborhoods  --  smallest: " + nbInfo[0] + "   largest: " + nbInfo[1]
						+ "   average: " + nbInfo[2] + "   h2-index: " + nbInfo[3]);
			}
		}
		// Level 3 statistics
		if (level > 2) {
			// Search for largest clique of graph
			if (!subCall) clique = CliqueAlgorithm.run(graph.getClone(), timeout, true);
			else {	// If the graph is a 2,t-club solution, we must map the IDs to 0..(size-1) for the copy
				/* We obtain a private copy of the graph with IDs mapped to ID range 0..(size-1) by calling
				 * getTwoNeighborhood for any vertex (2,t-clubs have diameter 2!). Granted, this method is a bit dirty. */
				RtcGraph copy = null;
				for (int anyVertex : graph.getVertices()) { copy = graph.getTwoNeighborhood(anyVertex); break; }
				clique = CliqueAlgorithm.run(copy, timeout, true);
			}
			maxClique = !CliqueAlgorithm.timedOut;
			if (maxClique) headerInfo.add("max. clique: "+clique.size());
				else headerInfo.add("max. clique: >= "+(clique.size()));
		}
		if (!silent) for (String line : headerInfo) System.out.println(line);
		header = headerInfo.toArray(new String[0]);
		return maxClique;
	}
	
	/** Analyzes given solution, i.e. 2,t-club. Outputs results to console, if 'silent' is not set true.
	 *  Additionally generates header information as string array for output files.
	 *  Returns the largest t for which the given graph is a solution.
	 * @param graph		The graph to be analyzed.
	 * @param level		Analysis level. 1 = Only basic information computable in linear time.
	 * 									2 = Additional information requiring higher order polynomial time.
	 * 									3 = Even more information requiring to solve NP-hard problems.
	 * @param t			Parameter t for which the graph is a maximum 2,t-club.
	 * @param silent	Whether output to the console shall be suppressed.
	 * @param header	Reference to string array in which header information shall be stored.
	 * @return			Maximum t for which solution is a 2,t-club.
	 * 					Header can be retrieved from public variable 'header' after call to this method.
	 */
	public static int analyzeSolution (RtcGraph graph, int level, int t, boolean silent, int timeout,
			boolean makeTable, Model usedModel) {
		// Get information from analyzeGraph method
		analyzeGraph(graph, level, "", silent, timeout, true);
		// Generate further information on solution
		ArrayList<String> furtherHeaderInfo = new ArrayList<String>();
		// Determine maximum t for which solution is a 2,t-club
		int bestT = graph.is2TClub(usedModel);
		int newV = -1;
		// Level 1 statistics
		if (compareSolution != null) {
			int commonV = getVertexSolutionDelta(graph);
			newV = graph.size() - commonV;
			int missingV = compareSolution.size() - commonV;
			furtherHeaderInfo.add("vertex solution delta  -- common: " + commonV + "  new: " + newV +
					"  missing: " + missingV);
		}	
		// Level 2 statistics
		int critVertices = 0;	// level 3 needs this value
		if (level > 1) {
			// Generate common neighborhood information
			int[] nbhInfo = getCommonNeighborhoodInfo(graph);
			furtherHeaderInfo.add("common neighborhoods  --  smallest: "+nbhInfo[0]+"  largest: "+nbhInfo[1]+
					"  average: "+nbhInfo[2]+"  cn-index: "+nbhInfo[3]);
			// Count critical edges, critical vertices and determine perfectness
			if (graph.size() < 1000 || level > 2) {		// gets "promoted" to level 3 for large solutions 
				int critEdges = countCriticalEdges(graph, bestT);
				critVertices = countCriticalVertices(graph, bestT);
				double perfectness = determinePerfectness(graph, bestT);
				furtherHeaderInfo.add((bestT>t?"For t="+bestT+":  ":"") + "critical vertices: "+critVertices+
						"   critical edges: "+critEdges+"   perfectness: "+format5.format(perfectness));
			}
		}
		// Level 3 statistics
		int degeneracy = t;
		if (level > 2) {
			// ...
			// Check degeneracy of solution
			Integer[] degenerateSet = null;
			for (degeneracy = bestT; degeneracy > 0; degeneracy--) {
				degenerateSet = checkIfDegenerate(graph, degeneracy, timeout);
				if (timedOut) break;
				if (degenerateSet != null) break;
			}
			if (timedOut) {
				furtherHeaderInfo.add("Could not determine if solution is "+degeneracy+"-degenerate, " +
						"due to algorithm timeout, ");
				furtherHeaderInfo.add("but could determine that it is not "+(degeneracy+1)+"-degenerate.");
			}
			else if (degeneracy > 0) {
				furtherHeaderInfo.add("");
				furtherHeaderInfo.add("This solution is "+degeneracy+"-degenerate, i.e. deleting only edges yields a");
				furtherHeaderInfo.add("complete bipartite graph with at least "+degeneracy+" "+ 
						(degeneracy==1?"vertex":"vertices") + " on either side.");
				if (degeneracy >= t) {
					furtherHeaderInfo.add("");
					if (degeneracy == 1) {
						furtherHeaderInfo.add("As degeneracy equals t = 1, this solution is strongly degenerate, meaning");
						furtherHeaderInfo.add("that it is just the closed 1-neighborhood of the maximum degree vertex");
						furtherHeaderInfo.add("of the original graph.");
					} else {
						furtherHeaderInfo.add("As degeneracy is at least t, this solution is strongly degenerate, meaning");
						furtherHeaderInfo.add("that it is just the largest common 1-neighborhood of any "+degeneracy+
								" vertices of the");
						furtherHeaderInfo.add("original graph plus the "+degeneracy+" vertices forming it.");
					}
				}
				furtherHeaderInfo.add("");
				String line = "Set of vertices causing "+degeneracy+"-degeneracy: ";
				for (int v : degenerateSet) {
					String newLine = line + " " + graph.getVertexName(v);
					if (newLine.length() > 80) {
						furtherHeaderInfo.add(line);
						newLine = "" + graph.getVertexName(v);
					}
					line = newLine;	
				}
				furtherHeaderInfo.add(line);
			}
			else if (degeneracy == 0) {
				furtherHeaderInfo.add("");
				furtherHeaderInfo.add("This solution is not degenerate, i.e. deleting only edges will not yield");
				furtherHeaderInfo.add("any complete bipartite graph, not even one with only one vertex on one side.");
			}
		}
		if (bestT > t) {
			furtherHeaderInfo.add(""); 
			furtherHeaderInfo.add("This solution is also a 2,t-club for t up to "+bestT+".");
		}
		if (!silent) for (String line : furtherHeaderInfo) System.out.println(line);
		// Combine header information
		String[] allHeaderInfo = new String[header.length + furtherHeaderInfo.size()];
		for (int i = 0; i < header.length; i++)
			allHeaderInfo[i] = header[i];
		for (int i = 0; i < furtherHeaderInfo.size(); i++)
			allHeaderInfo[i+header.length] = furtherHeaderInfo.get(i);
		header = allHeaderInfo;
		// Pass information to LaTeX TableMaker, if output = true and level = 3
		if (makeTable && level == 3) {
			long pairs = (long)(graph.size() - 1) * graph.size() / 2;
			double density = (double)(graph.getEdgeCount()) / pairs;
			TableMaker.addSolutionInfo(t, graph.size(), graph.getEdgeCount(), density, newV, degeneracy);
		}	
		return bestT;
	}
	
	public static void setCompareSolution (RtcGraph solution) {
		compareSolution = solution;
	}
	
	// Returns number of common vertices in new and old solution
	private static int getVertexSolutionDelta (RtcGraph solution) {
		HashSet<String> oldVertices = compareSolution.getVertexNames();
		HashSet<String> newVertices = solution.getVertexNames();
		int common = 0;
		for (String n : oldVertices)
			if (newVertices.contains(n))
				common++;
		return common;
	}
	
	/**
	 * 
	 * @param graph
	 * @return Integer array of four components: 	1. minimum degree
	 * 												2. maximum degree
	 * 												3. average degree (rounded down)
	 * 												4. h-index
	 */
	private static int[] getDegreeInfo (RtcGraph graph) {
		int[] solutions = new int[4];
		solutions[0] = -1;
		int degreeSum = 0;
		int degreeCount[] = new int[graph.size()];
		for (int v : graph.getVertices()) {
			if (solutions[0] == -1 || graph.degree(v) < solutions[0])
				solutions[0] = graph.degree(v);
			if (graph.degree(v) > solutions[1])
				solutions[1] = graph.degree(v);
			degreeSum += graph.degree(v);
			degreeCount[graph.degree(v)]++;
		}
		solutions[2] = degreeSum / graph.size();
		if (degreeSum / 1.0 / graph.size() % 1 >= 0.5) solutions[2]++;
		int degreePtr = graph.size();
		int verticesCounted = 0;
		while (degreePtr > verticesCounted) {
			degreePtr--;
			verticesCounted += degreeCount[degreePtr];
		}
		solutions[3] = degreePtr;
		return solutions;
	}
	
	/**
	 * 
	 * @param graph
	 * @return Integer array of four components: 	1. smallest 2-neighborhood
	 * 												2. largest 2-neighborhood
	 * 												3. average 2-neighborhood
	 * 												4. h2-index
	 */
	private static int[] get2NeighborhoodInfo (RtcGraph graph) {
		int[] solutions = new int[4];
		solutions[0] = -1;
		int nbSum = 0;
		int nbCount[] = new int[graph.size()];
		for (int v : graph.getVertices()) {
			int nb = graph.sizeOfTwoNeighborhood(v, false);
			if (solutions[0] == -1 || nb < solutions[0])
				solutions[0] = nb;
			if (nb > solutions[1])
				solutions[1] = nb;
			nbSum += nb;
			nbCount[nb]++;
		}
		solutions[2] = nbSum / graph.size();
		if (nbSum / 1.0 / graph.size() % 1 >= 0.5) solutions[2]++;
		int nbPtr = graph.size();
		int verticesCounted = 0;
		while (nbPtr > verticesCounted) {
			nbPtr--;
			verticesCounted += nbCount[nbPtr];
		}
		solutions[3] = nbPtr;
		return solutions;
	}
	
	/**
	 * 
	 * @param graph
	 * @return Integer array of four components: 	1. smallest common neighborhood
	 * 												2. largest common neighborhood
	 * 												3. average common neighborhood
	 * 												4. CN-index (= at least that many vertices have that many common
	 * 															 neighbors with all non-adjacent vertices)
	 */
	private static int[] getCommonNeighborhoodInfo (RtcGraph graph) {
		int[] solutions = new int[4];
		solutions[0] = -1;
		int nbSum = 0;
		int nbCount[] = new int[graph.size()];
		for (int v : graph.getVertices()) {
			int minNB = graph.size() - 1;		// Vertices adjacent to all others keep this value
			for (int w : graph.getVertices()) if (w > v) {
				int nb = graph.countCommonNeighbors(v, w);
				if (!graph.adjacent(v, w) && nb < minNB) minNB = nb;
				if (w > v) {
					if (solutions[0] == -1 || nb < solutions[0])
						solutions[0] = nb;
					if (nb > solutions[1])
						solutions[1] = nb;
					nbSum += nb;
				}
			}
			nbCount[minNB]++;
		}
		solutions[2] = nbSum / (graph.size() * (graph.size() - 1) / 2);
		if (nbSum / 1.0 / (graph.size() * (graph.size() - 1) / 2) % 1 >= 0.5) solutions[2]++;
		int nbPtr = graph.size();
		int verticesCounted = 0;
		while (nbPtr > verticesCounted) {
			nbPtr--;
			verticesCounted += nbCount[nbPtr];
		}
		solutions[3] = nbPtr;
		return solutions;
	}
	
	private static int countCriticalEdges (RtcGraph graph, int t) {
		int critical = 0;
		for (int v : graph.getVertices())
			for (int w : graph.getNeighbors(v)) if (w > v)
				if (graph.countCommonNeighbors(v, w) < t)
					critical++;
		return critical;
	}
	
	private static int countCriticalVertices (RtcGraph graph, int t) {
		HashSet<Integer> critical = new HashSet<Integer>();
		for (int v : graph.getVertices()) for (int w : graph.getVertices())
			if (w > v && !graph.adjacent(v, w))
				if (graph.countCommonNeighbors(v, w) == t)
					critical.addAll(graph.getCommonNeighbors(v, w));
		return critical.size();
	}
	
	private static double determinePerfectness (RtcGraph graph, int t) {
		long neighbors = 0;
		long perfectAmount = 0;
		for (int v : graph.getVertices()) for (int w : graph.getVertices()) if (w > v) {
			neighbors += graph.countCommonNeighbors(v, w);
			if (!graph.adjacent(v, w)) perfectAmount += t;
		}
		return (double)(perfectAmount)/neighbors;
	}
	
	/** A solution is d-degenerate, if n-d vertices are common neighbors of d vertices.
	 *  @return Null, if not having demanded degeneracy. Array of vertices causing strong degeneracy otherwise. */
	private static Integer[] checkIfDegenerate (RtcGraph graph, int degeneracy, int timeout) {
		timedOut = false;
		start = System.currentTimeMillis();
		int highDegree = graph.size() - degeneracy;
		if (highDegree < degeneracy) return null;
		HashSet<Integer> candidates = new HashSet<Integer>();
		for (int v : graph.getVertices())
			if (graph.degree(v) >= highDegree) {
				if (degeneracy == 1) {
					Integer[] array = new Integer[1];
					array[0] = v;
					return array;
				} else candidates.add(v);
			}
		if (candidates.size() < degeneracy) return null;
		// Build base set, i.e. set of vertices seeing all vertices
		HashSet<Integer> baseSet = new HashSet<Integer>();
		for (int v : candidates) {
			if (graph.degree(v) == graph.size() - 1) baseSet.add(v);
			if (baseSet.size() == degeneracy) return baseSet.toArray(new Integer[0]);
		}
		candidates.removeAll(baseSet);
		Integer[] degenerateSet = null;
		for (int c : candidates) {
			Integer[] group = new Integer[1+baseSet.size()];
			group[0] = c;
			int ptr = 1;
			for (int bv : baseSet) group[ptr++] = bv;
			HashSet<Integer> neighbors = new HashSet<Integer>();
			for (int n : graph.getNeighbors(c)) neighbors.add(n);
			HashSet<Integer> remainingCandidates = new HashSet<Integer>();
			for (int rc : candidates) if (rc != c) remainingCandidates.add(rc);
			degenerateSet = checkGroup(graph, degeneracy, highDegree, group, neighbors, remainingCandidates, timeout);
			if (degenerateSet != null) break;
			if (timedOut) return null;
		}
		return degenerateSet;
	}
	
	private static Integer[] checkGroup (RtcGraph graph, int degeneracy, int bound, Integer[] group, 
			HashSet<Integer> neighbors, HashSet<Integer> cand, int timeout) {
		HashSet<Integer[]> newGroups = new HashSet<Integer[]>();
		HashMap<Integer[], Integer> largest = new HashMap<Integer[], Integer>();
		HashMap<Integer[], HashSet<Integer>> neighborMap = new HashMap<Integer[], HashSet<Integer>>();
		
		// Expand group by adding all candidates not seen by all vertices in group
		HashSet<Integer> expansionSet = new HashSet<Integer>();
		for (int c : cand) if (!neighbors.contains(c)) expansionSet.add(c);
		if (!expansionSet.isEmpty()) {
			if (expansionSet.size() + group.length > degeneracy) return null;
			Integer[] expandedGroup = new Integer[group.length + expansionSet.size()];
			for (int ptr = 0; ptr < group.length; ptr++)
				expandedGroup[ptr] = group[ptr];
			int ptr = group.length;
			for (int e : expansionSet) expandedGroup[ptr++] = e;
			HashSet<Integer> remainingNeighbors = new HashSet<Integer>();
			for (int n : neighbors) {
				boolean adjacentToAll = true;
				for (int e : expansionSet)
					if (!graph.adjacent(n, e)) {adjacentToAll = false; break;}
				if (adjacentToAll) remainingNeighbors.add(n);
			}
			if (remainingNeighbors.size() < bound) return null;
			if (expandedGroup.length == degeneracy) return expandedGroup;
			HashSet<Integer> remainingCandidates = new HashSet<Integer>();
			for (int c : cand) if (!expansionSet.contains(c)) remainingCandidates.add(c);
			return checkGroup(graph, degeneracy, bound, expandedGroup, remainingNeighbors, remainingCandidates, timeout);
		}
		
		long timePassed = System.currentTimeMillis() - start;
		if (timeout != -1 && timePassed > timeout * 1000) {
			timedOut = true; return null;
		}
		
		// Build base set, i.e. set of vertices seeing all vertices in set 'neighbors'
		HashSet<Integer> baseSet = new HashSet<Integer>();
		HashSet<Integer> candidates = new HashSet<Integer>();	// need copy because other branches still need original
		for (int v : cand) {
			boolean seesThemAll = true;
			for (int n : neighbors)
				if (!graph.adjacent(v, n)) {seesThemAll = false; break;}
			if (seesThemAll) baseSet.add(v);
			else candidates.add(v);
			if (baseSet.size() + group.length == degeneracy) {
				for (int g : group) baseSet.add(g);
				return baseSet.toArray(new Integer[0]);
			}
		}
		int oldLen = group.length;
		int len = group.length+1+baseSet.size();
			for (int c : candidates) {
				// Only continue, if c is larger than all members of the old group.
				// This avoids creating any subset twice.
				boolean cIsLargest = true;
				for (int g : group) if (g >= c) cIsLargest = false;
				if (!cIsLargest) continue;
				Integer[] newGroup = new Integer[len];
				for (int i = 0; i < oldLen; i++)
					newGroup[i] = group[i];
				int ptr = oldLen;
				for (int bv : baseSet) newGroup[ptr++] = bv;
				newGroup[len-1] = c;
				HashSet<Integer> newNeighbors = new HashSet<Integer>();
				for (int n : neighbors) if (graph.adjacent(c, n))
					newNeighbors.add(n);
				if (newNeighbors.size() >= bound) {
					if (len >= degeneracy) return newGroup;
					newGroups.add(newGroup);
					largest.put(newGroup, c);
					neighborMap.put(newGroup, newNeighbors);
				}	
			}
			if (len == degeneracy) return null;
			Integer[] degenerateSet = null;
			for (Integer[] g : newGroups) {
				int max = largest.get(g);
				HashSet<Integer> remainingCandidates = new HashSet<Integer>();
				for (int c : candidates)
					if (c > max) remainingCandidates.add(c);
				if (remainingCandidates.size() >= degeneracy - g.length)
					degenerateSet = checkGroup(graph, degeneracy, bound, g, neighborMap.get(g), remainingCandidates, timeout);
				if (degenerateSet != null) return degenerateSet;
				if (timedOut) break;
			}	
			return null;
		}
	
	/**	Combines input graph header and solution graph header into one header for the output file.
	 * 
	 * @param graphHeader		The header information for the input graph.
	 * @param solutionHeader	The header information for the output graph.
	 * @param t					The parameter t for which the output graph is a solution.
	 * @param graphName			The name of the input graph, usually the file name.
	 * @param optimal			Whether the solution is optimal.
	 * @param timeout			Timeout of the algorithm in seconds.
	 * @return					Complete header for the output file.
	 */
	public static String[] combineHeaders (String[] graphHeader, String[] solutionHeader, int t, String graphName,
			boolean optimal, int timeout) {
		int opt = optimal?0:2;
		String[] header = new String[2 + opt];
		if (graphHeader != null && solutionHeader != null)
			header = new String[graphHeader.length + solutionHeader.length + 4 + opt];
		else if (graphHeader != null)
			header = new String[graphHeader.length + 3 + opt];
		else if (solutionHeader != null)
			header = new String[solutionHeader.length + 3 + opt];
		if (optimal)
			header[0] = "# This graph is a maximum size 2,"+t+"-club of the graph \""+graphName+"\".";
		else {
			header[0] = "# This graph is the largest found 2,"+t+"-club of the graph \""+graphName+"\".";
			if (timeout < 60)
				header[1] = "# The 2,t-club algorithm timed out after " + timeout + " seconds,";
			else if (timeout < 3600)
				header[1] = "# The 2,t-club algorithm timed out after " + timeout/60 + " minutes and "+
							timeout%60 + " seconds,";
			else if (timeout < 24 * 3600)
				header[1] = "# The 2,t-club algorithm timed out after " + timeout/3600 + " hours, "+
							(timeout/60)%60 + " minutes and " + timeout%60 + " seconds,";
			else header[1] = "# The 2,t-club algorithm timed out after " + timeout/(24*3600) + " days," + 
							((timeout/3600)%24) + " hours and " + (timeout/60)%60 + " minutes,";		
			header[2] = "# therefore this solution is not guaranteed to be optimal.";
		}
		if (graphHeader != null) {
			header[1+opt] = "# ";
			for (int line = 0; line < graphHeader.length; line++)
				header[2+line+opt] = "# " + graphHeader[line];
		}
		if (solutionHeader != null) {
			int graphHeaderLength = (graphHeader == null)?0:(graphHeader.length + 1);
			header[1+graphHeaderLength+opt] = "# ";
			for (int line = 0; line < solutionHeader.length; line++)
				header[2+graphHeaderLength+line+opt] = "# " + solutionHeader[line];
		}
		header[header.length - 1] = "# ";
		return header;
	}
	
}
