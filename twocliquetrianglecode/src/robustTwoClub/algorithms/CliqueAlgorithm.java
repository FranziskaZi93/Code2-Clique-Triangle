package robustTwoClub.algorithms;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.ListIterator;
import java.util.Stack;

import robustTwoClub.graph.RtcGraph;


public class CliqueAlgorithm {

	private static HashSet<Integer> solution;
	private static int lowerBound;
	
	private static int currentKernel;
	private static int kernelSize;
	private static HashMap<Integer, HashSet<Integer>> adjacency;
	private static int[] mapBack;
	private static HashMap<Integer, Integer> map;
	private static HashSet<Integer> fixed;
	
	private static Stack<Integer> rollbackID;
	private static Stack<HashSet<Integer>> rollbackAdj;
	
	private static int timeout;
	private static long start;
	private static boolean timeRunOut;
	
	private static boolean silent;
	
	public static boolean timedOut;
	
	private static final DecimalFormat tf = new DecimalFormat("0.000");
	
	/**
	 * 
	 * @param graph		Graph on which to run algorithm. Algorithm reduces graph, so pass copy!
	 * 					Vertex IDs are required to be in range 0..(graph size - 1)!
	 * @param timeout	Allowed running time in seconds. Pass -1 for no limit.
	 * @return
	 * @throws GraphException
	 */
	public static HashSet<Integer> run (RtcGraph graph, int timeOut, boolean beSilent) {
		timeout = timeOut;
		timeRunOut = false;
		start = System.nanoTime();
		silent = beSilent;
		ArrayList<Integer> toSort = new ArrayList<Integer>();
		int[] kernelSizes = new int[graph.size()];
		for (int v : graph.getVertices()) {
			toSort.add(v);
			kernelSizes[v] = graph.degree(v);
		}	
		ArrayList<Integer> kernelOrder = kernelQuicksort(toSort, kernelSizes);
		ListIterator<Integer> li = kernelOrder.listIterator();
		lowerBound = 0; solution = null;
		while (li.hasNext()) {
			if (timeout != -1 && System.nanoTime() - start > timeout * 1000000000) {
				timeRunOut = true; break;
			}
			int next = li.next();
			if (graph.degree(next) >= lowerBound)
				searchOnKernel(next, graph);
			graph.deleteVertex(next);
		}
		long time = System.nanoTime() - start;
		if (!silent) {
			if (timeRunOut)
				System.out.println("\rClique algorithm reached timeout of "+timeout+ " seconds.                              ");
			else System.out.println("\rClique algorithm running time: "+tf.format(time/1000000000.0)+ " seconds              ");
			if (solution == null)
				System.out.println("Found no clique due to timeout.");
			else if (timeRunOut)
				System.out.println("Found clique of size "+solution.size()+", which is NOT GUARANTEED to be"+
						" maximum due to timeout.");
			else System.out.println("Found maximum clique of size "+solution.size());
			HashSet<String> names = new HashSet<String>();
			if (solution != null) {
				for (int v: solution)
					names.add(graph.getVertexName(v));
				System.out.println("Vertices in the clique: "+solution.toString());
			}	
		}
		timedOut = timeRunOut;
		return solution;
	}
	
	private static ArrayList<Integer> kernelQuicksort (ArrayList<Integer> toSort, int[] kernelSizes) {
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
		if (!smallerList.isEmpty()) solution.addAll(kernelQuicksort(smallerList, kernelSizes));
		solution.addAll(equalsList);
		if (!greaterList.isEmpty()) solution.addAll(kernelQuicksort(greaterList, kernelSizes));
		return solution;
	}
	
	private static void searchOnKernel (int kernel, RtcGraph graph) {
		currentKernel = kernel;
		mapBack = new int[graph.degree(kernel)];
		map = new HashMap<Integer, Integer>();
		adjacency = new HashMap<Integer, HashSet<Integer>>();
		fixed = new HashSet<Integer>();
		rollbackID = new Stack<Integer>();
		rollbackAdj = new Stack<HashSet<Integer>>();
		int counter = 0;
		for (int v : graph.getNeighbors(kernel)) { 
				mapBack[counter] = v;
				map.put(v, counter);
				adjacency.put(counter, new HashSet<Integer>());
				counter++;
			}
		for (int i = 0; i < adjacency.size(); i++)
			for (int v : graph.getNeighbors(mapBack[i]))
				if (graph.adjacent(v, kernel))
					adjacency.get(i).add(map.get(v));
		kernelSize = adjacency.size();
		// Purge vertices with less than lowerBound minus one neighbors in local adjacency structure
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int i = 0; i < adjacency.size(); i++)
			if (adjacency.get(i).size() < lowerBound - 1)
				deleteSet.add(i);
		while (!deleteSet.isEmpty()) {
			int next = 0;								// dummy initialization
			for (int v : deleteSet) {next = v; break;}	// grab any vertex in delete set
			deleteSet.remove(next);
			for (int nb : adjacency.get(next)) {
				adjacency.get(nb).remove(next);
				if (adjacency.get(nb).size() == lowerBound - 2)
					deleteSet.add(nb);
			}
			adjacency.remove(next);
			if (adjacency.size() - deleteSet.size() < lowerBound) return;
		}
		// Start branching
		branch();
		return;
	}
	
	private static void branch () {
		
		if (timeout != -1 && System.nanoTime() - start > timeout * 1000000000) timeRunOut = true;
		if (timeRunOut) return;
		
		// Determine minimum degree vertex
		int minDegVertex = -1;
		for (int v : adjacency.keySet())
			if (!fixed.contains(v))
				if (minDegVertex == -1 || adjacency.get(v).size() < adjacency.get(minDegVertex).size())
					minDegVertex = v;
		
		// Check for solution
		if (minDegVertex == -1 || adjacency.get(minDegVertex).size() == adjacency.size() - 1) {
			lowerBound = adjacency.size() + 1;
			solution = new HashSet<Integer>();
			for (int v : adjacency.keySet())
				solution.add(mapBack[v]);
			solution.add(currentKernel);
			if (!silent) System.out.print("\rCurrently best solution size: "+solution.size()+"  ");
			return;
		}
		
		// Check color bound
		if (fixed.size() + colorBound() < lowerBound) return;		
		
		// First branch over deleting minimum degree vertex then over fixing it
		int rollbackPoint = -1;
		if (!rollbackID.isEmpty()) rollbackPoint = rollbackID.peek();
		boolean branch = delete(minDegVertex);
		if (branch) branch();
		rollback(rollbackPoint);
		rollbackPoint = -1;
		if (!rollbackID.isEmpty()) rollbackPoint = rollbackID.peek();
		branch = fix(minDegVertex);
		if (branch) branch();
		rollback(rollbackPoint);
		fixed.remove(minDegVertex);
		return;
	}
	
	private static boolean fix (int vertex) {
		fixed.add(vertex);
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int v : adjacency.keySet())
			if (v != vertex && !adjacency.get(vertex).contains(v))
				deleteSet.add(v);
		return delete(deleteSet);
	}
	
	private static boolean delete (int vertex) {
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		deleteSet.add(vertex);
		return delete(deleteSet);
	}
	
	private static boolean delete (HashSet<Integer> deleteSet) {
		while (!deleteSet.isEmpty()) {
			int next = 0;								// dummy initialization
			for (int v : deleteSet) {next = v; break;}	// grab any vertex in delete set
			rollbackID.push(next);
			HashSet<Integer> adjEntry = adjacency.remove(next);
			rollbackAdj.push(adjEntry);
			for (int v : adjEntry) {
				adjacency.get(v).remove(next);
				if (adjacency.get(v).size() == lowerBound - 2)
					deleteSet.add(v);
			}
			deleteSet.remove(next);
			if (adjacency.size() - deleteSet.size() < lowerBound) return false;
		}
		return true;
	}
	
	
	private static void rollback (int rollbackPoint) {
		while (!rollbackID.isEmpty() && rollbackID.peek() != rollbackPoint) {
			int id = rollbackID.pop();
			HashSet<Integer> adjEntry = rollbackAdj.pop();
			adjacency.put(id, adjEntry);
			for (int v : adjEntry)
				adjacency.get(v).add(id);
		}
	}
	
	private static int colorBound () {
		int colors[] = new int[kernelSize+1];
		int highestColor = 0;
		for (int v : adjacency.keySet()) 
			if (!fixed.contains(v)) {
				boolean colorsSeen[] = new boolean[kernelSize+1];
				for (int w : adjacency.get(v))
					if (colors[w] != 0)
						colorsSeen[colors[w]] = true;
				int color = 1;
				while (colorsSeen[color])
					color++;
				colors[v] = color;
				if (color > highestColor)
					highestColor = color;
			}
		return highestColor;	// upper bound on clique size
	}
	
}
