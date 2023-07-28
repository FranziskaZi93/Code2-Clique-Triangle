package robustTwoClub.experiments;

import java.util.ArrayList;
import java.util.HashMap;

import robustTwoClub.algorithms.RobustTwoClubAlgorithm.Model;
import robustTwoClub.graph.RtcGraph;

public class CNBoundedGenerator {
	
	private static RtcGraph graph;
	private static int[][] cnb;
	private static HashMap<Integer, ArrayList<Integer>> conflicts;
	private static int conflictNr;
	
	public static void main (String[] args) {
		bruteForceSearch(16, 2, 2, true, false, true, 100);
		//System.out.println(experimentallyDeriveGraphSizeLowerBound (2, 16, false, 1, true));
	}
	
	/** Derives lower bound for maximum size of a 2,t-club containing at most b common vertices per vertex pair. */
	public static int experimentallyDeriveGraphSizeLowerBound (int t, int b, boolean tFree, int limit, boolean outputBest) {
		int currentSize = 1;
		int largestKnownPossibleSize = 1;
		int smallestKnownDifficultSize = -1;
		RtcGraph bestGraph = null;
		while (smallestKnownDifficultSize == -1) {
			currentSize *= 2;
			if (currentSize < t) currentSize = t;
			RtcGraph graph = bruteForceSearch(currentSize, t, b, tFree, true, false, limit);
			boolean foundOne = (graph != null);
			if (foundOne && (bestGraph == null || (graph.size() > bestGraph.size()))) bestGraph = graph;
			if (foundOne) largestKnownPossibleSize = currentSize;
			else smallestKnownDifficultSize = currentSize;
		}
		if (smallestKnownDifficultSize <= 36) smallestKnownDifficultSize = 36; // around small values there can be gaps!
		while (smallestKnownDifficultSize - largestKnownPossibleSize > 1) {
			currentSize = (largestKnownPossibleSize + smallestKnownDifficultSize) / 2;
			if (smallestKnownDifficultSize <= 36) currentSize = smallestKnownDifficultSize - 1; 
			RtcGraph graph = bruteForceSearch(currentSize, t, b, tFree, true, false, limit);
			boolean foundOne = (graph != null);
			if (foundOne && (bestGraph == null || (graph.size() > bestGraph.size()))) bestGraph = graph;
			if (foundOne) largestKnownPossibleSize = currentSize;
			else smallestKnownDifficultSize = currentSize;
		}
		if (outputBest) graph.writeToFile("output", "best"+b+"CNbounded"+(tFree?"T-Free":"")+"2,"+t+"Club", null);
		return largestKnownPossibleSize;
	}
	
	/** Returns a graph with specified properties, if one has been found. */
	public static RtcGraph bruteForceSearch (int n, int t, int b, boolean tFree, boolean silent, boolean out, int limit) {
		RtcGraph graph = generateCNBounded2TClub(n, t, b, tFree);
		int runs = 1;
		while (graph == null && runs < limit) {
			if (!silent) System.out.println("Failed to construct "+b+"-CN-bounded "+(tFree?"triangle-free ":"")+
					"2,"+t+"-club of "+n+" vertices.");
			graph = generateCNBounded2TClub(n, t, b, tFree);
			runs++;
		}
		if (graph == null) return null;
		if (graph.size() > t && graph.is2TClub(Model.HEREDITARY) < t) 
			System.out.println("WARNING: Generated graph is no 2,"+t+"-club!");
		boolean cnBounded = true;
		for (int v : graph.getVertices())
			for (int w: graph.getVertices())
				if (v != w)
					if (graph.countCommonNeighbors(v, w) > b)
						cnBounded = false;
		if (!cnBounded)	System.out.println("WARNING: Generated graph is not "+b+"-CN-bounded!");
		boolean triangleFree = true;
		if (tFree)
			for (int v : graph.getVertices())
				for (int w : graph.getNeighbors(v))
					if (graph.countCommonNeighbors(v, w) > 0)
						triangleFree = false;
		if (!triangleFree)	System.out.println("WARNING: Generated graph is not triangle-free!");
		if (!silent) System.out.println("Generation complete!");
		if (out) graph.writeToFile("output", ""+b+"CNboundedT"+t+"N"+n+(tFree?"T-Free":""), null);
		return graph;
	}


	public static RtcGraph generateCNBounded2TClub (int size, int t, int b, boolean tFree) {
		graph = new RtcGraph();
		for (int i = 0; i < size; i++)
			graph.addVertex(i, ""+i);
		cnb = new int[size][size];
		conflicts = new HashMap<Integer, ArrayList<Integer>>();
		for (int v = 0; v < size; v++) {
			conflicts.put(v, new ArrayList<Integer>());
			for (int w = 0; w < size; w++)
				if (w != v)
					conflicts.get(v).add(w);
		}
		conflictNr = (size * (size - 1)) / 2;
		while (conflictNr != 0) {
			// Pick random vertex with remaining conflicts
			int v = (int)(Math.random() * size);
			while (conflicts.get(v).size() == 0)
				v = (int)(Math.random() * size);
			// Pick random vertex from conflicts of v
			int wPos = (int)(Math.random() * conflicts.get(v).size());
			int w = conflicts.get(v).get(wPos);
			// Try to add an edge for v and w
			if (tryToAddEdge(v, w, t, b, tFree)) continue;
			// Else try to add new common neighbors
			for (int c : graph.getVertices()) {
				if (c == v || c == w) continue;
				if (graph.adjacent(v, c) && graph.adjacent(w, c)) continue;
				if (graph.adjacent(v, c)) tryToAddEdge(w, c, t, b, tFree);
				else if (graph.adjacent(w, c)) tryToAddEdge(v, c, t, b, tFree);
				else if (tryToAddCommonNeighbor(v, w, c, b, tFree)) {
					tryToAddEdge(v, c, t, b, tFree);
					tryToAddEdge(w, c, t, b, tFree);
				}
				if (cnb[v][w] == t) break;
			}
			if (cnb[v][w] < t) return null;
		}
		return graph;
	}
	
	private static boolean tryToAddEdge (int v, int w, int t, int b, boolean tFree) {
		// Check whether we may add an edge for v and w
		boolean mayAddEdge = true;
		if (graph.countCommonNeighbors(v, w) > 0 && tFree)
			mayAddEdge = false;
		for (int x : graph.getNeighbors(v))
			if (cnb[w][x] == b)
				mayAddEdge = false;
		for (int x : graph.getNeighbors(w))
			if (cnb[v][x] == b)
				mayAddEdge = false;
		// Add edge and update information if we may
		if (mayAddEdge) {
			if (cnb[v][w] < t) conflictNr--;
			for (int x : graph.getNeighbors(v)) {
				cnb[w][x]++;
				cnb[x][w]++;
				if (cnb[w][x] == t && !graph.adjacent(w, x)) {
					conflictNr--;
					conflicts.get(w).remove(new Integer(x));
					conflicts.get(x).remove(new Integer(w));
				}
			}
			for (int x : graph.getNeighbors(w)) {
				cnb[v][x]++;
				cnb[x][v]++;
				if (cnb[v][x] == t && !graph.adjacent(v, x)) {
					conflictNr--;
					conflicts.get(v).remove(new Integer(x));
					conflicts.get(x).remove(new Integer(v));
				}
			}
			conflicts.get(v).remove(new Integer(w));
			conflicts.get(w).remove(new Integer(v));
			graph.addEdge(v, w);		
		}
		return mayAddEdge;
	}
	
	private static boolean tryToAddCommonNeighbor (int v, int w, int c, int b, boolean tFree) {
		// Check whether we may add both edges (v,c) and (w,c)
		if (tFree) {
			if (graph.adjacent(v, w)) return false;
			if (graph.countCommonNeighbors(v, c) > 0) return false;
			if (graph.countCommonNeighbors(w, c) > 0) return false;
		}
		for (int x : graph.getNeighbors(v))
			if (cnb[c][x] == b) return false;
		for (int x : graph.getNeighbors(w))
			if (cnb[c][x] == b) return false;
		for (int x : graph.getNeighbors(c)) {
			if (cnb[v][x] == b) return false;
			if (cnb[w][x] == b) return false;
		}
		for (int x : graph.getNeighbors(v))
			if (graph.getNeighbors(w).contains(x))
				if (cnb[c][x] == b-1) return false;
		return true;
	}
	
}
