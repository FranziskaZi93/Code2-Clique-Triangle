package robustTwoClub.graph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.*;

import robustTwoClub.algorithms.RobustTwoClubAlgorithm.Model;

public class RtcGraph {

	public static long	flowConstructGraph, flowComputePath;

	private HashMap<Integer, String> idMap;					// remembers vertex names for outputting solution graphs
	private HashSet<Integer> nodes;							// set of vertices of the graph
	private HashMap<Integer, HashSet<Integer>> adjacency;	// sparse adjacency matrix of the graph

	/** Loads a graph from a DIMACS or METIS file.
	 *
	 * @param filename Name of the file to be read.
	 * @param type 0 = DIMACS, 1 = METIS, 2 = DIMACS clq
	 * @throws GraphException
	 */
	public RtcGraph (String filename, int type) {

		idMap = new HashMap<Integer, String>();
		nodes = new HashSet<Integer>();
		adjacency = new HashMap<Integer, HashSet<Integer>>();
		buildGraphFromFile(filename, type);
	}

	/**	Tries to guess whether input file has DIMACS or METIS format.
	 *  Guess will be nonsense, if file is neither.
	 *
	 * @param filename	Name of the input file for which the type shall be guessed.
	 * @return			Whether the input seems to be DIMACS. (METIS otherwise)
	 */
	public static int guessInputType (String filename) {
		if (filename.endsWith(".dimacs")) return 0;	// I guess it is a DIMACS graph
		if (filename.endsWith(".edges")) return 0;	// I guess it is a DIMACS graph
		if (filename.endsWith(".mtx")) return 0;
		if (filename.endsWith(".graph")) return 1;	// I guess it is a METIS graph
		if (filename.endsWith(".clq")) return 2;	// I guess it is a DIMACS clq graph
		if (filename.endsWith(".col")) return 2;
		System.out.println("Cannot determine type of input file.");
		System.out.println("If file is either a valid DIMACS or Metis file, please specify file type");
		System.out.println("by appending options '-dimacs' or '-metis' to call of program.");
		System.exit(0); // Abort program if input type can't be guessed from extension
		return 0;	// Dummy return
	}

	public int size () {
		return nodes.size();
	}

	/** Returns whether the given vertex ID belongs to the graph. */
	public boolean contains (int v) {
		return nodes.contains(v);
	}

	public int degree (int v) {
		return adjacency.get(v).size();
	}

	/** Returns whether vertices v and w are adjacent. */
	public boolean adjacent (int v, int w) {
		return adjacency.get(v).contains(w);
	}

	public HashSet<Integer> getVertices () {
		return nodes;
	}

	public int getEdgeCount () {
		int edges = 0;
		for (int v : nodes)
			edges += adjacency.get(v).size();
		edges /= 2;
		return edges;
	}

	public String edgesToString ()
	{
		return adjacency.toString();
	}

	public HashSet<Integer> getNeighbors (int v) {
		return adjacency.get(v);
	}

	public HashSet<Integer> getTwoNeighbors (int v) {
		HashSet<Integer> nb = new HashSet<Integer>();
		for (int w : adjacency.get(v))
			nb.addAll(adjacency.get(w));
		nb.addAll(adjacency.get(v));
		nb.remove(v);
		return nb;
	}

	public int countCommonNeighbors (int v, int w) {
		int neighbors = 0;
		for (int x : adjacency.get(v))
			if (adjacency.get(w).contains(x))
				neighbors++;
		return neighbors;
	}

	public int countCommonNeighbors (Integer[] vertices) {
		int neighbors = 0;
		for (int x : adjacency.get(vertices[0])) {
			boolean compatible = true;
			for (int i = 1; i < vertices.length; i++)
				if (!adjacency.get(vertices[i]).contains(x))
				{compatible = false; break;}
			if (compatible) neighbors++;
		}
		return neighbors;
	}

	public HashSet<Integer> getCommonNeighbors (int v, int w){
		if(adjacency.get(w).size() < adjacency.get(v).size()){
			int x = v;
			v = w;
			w = x;
		}
		HashSet<Integer> neighbors = new HashSet<Integer>();
		for (int x : adjacency.get(v))
			if (adjacency.get(w).contains(x))
				neighbors.add(x);
		return neighbors;
	}

	public int sizeOfTwoNeighborhood (int x, boolean countX) {
		HashSet<Integer> vertices = new HashSet<Integer>();
		for (int v : adjacency.get(x)) {
			vertices.add(v);
			for (int w : adjacency.get(v))
				vertices.add(w);
		}
		if (!countX) vertices.remove(x);
		return vertices.size();
	}

	/** Creates a subgraph of given vertices and remaps vertex IDs to range 0 .. (subgraph size - 1) */
	public RtcGraph getSubgraph (HashSet<Integer> vertices, int center) {
		RtcGraph subgraph = new RtcGraph();
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapBack = new HashMap<Integer, Integer>();
		// Add center
		map.put(center, 0);
		map.put(0, center);
		subgraph.addVertex(0);
		int id = 1;
		// Add vertices
		for (int v : vertices) if (v != center) {
			map.put(v, id);
			mapBack.put(id, v);
			subgraph.addVertex(id++);
		}
		// Add edges
		for (int v : vertices) for (int w : adjacency.get(v))
			if (w > v && vertices.contains(w))
				subgraph.addEdge(map.get(v), map.get(w));
		// Construct new idMap
		HashMap<Integer, String> newIdMap = new HashMap<Integer, String>();
		for (int i = 0; i < subgraph.size(); i++)
			newIdMap.put(i, idMap.get(mapBack.get(i)));
		subgraph.idMap = newIdMap;
		return subgraph;
	}

	/** Returns the 2-neighborhood of some vertex as a new graph.
	 *
	 * @param c Vertex for which to obtain the 2-neighborhood.
	 * @return The new graph.
	 */
	public RtcGraph getTwoNeighborhood (Integer c) {
		RtcGraph g = new RtcGraph();
		// We map the original integer IDs to new ones in order to minimize the range of vertex IDs
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapBack = new HashMap<Integer, Integer>();
		g.addVertex(0);										// add center vertex
		map.put(c, 0);
		mapBack.put(0, c);
		int index = 1;
		for (int v : adjacency.get(c)) {
			if (!map.containsKey(v)) {
				g.addVertex(index);
				map.put(v, index);
				mapBack.put(index, v);
				index++;
			}
			for (int w : adjacency.get(v)) {
				if (!map.containsKey(w)) {
					g.addVertex(index);
					map.put(w, index);
					mapBack.put(index, w);
					index++;
				}
			}
		}
		for (int v : g.nodes)
			for (int w : adjacency.get(mapBack.get(v)))
				if (map.containsKey(w) && map.get(w) > v)
					g.addEdge(v, map.get(w));

		// Construct idMap for new graph, mapping new integer IDs to original string identifiers
		HashMap<Integer, String> newIdMap = new HashMap<Integer, String>();
		for (int i = 0; i < g.size(); i++)
			newIdMap.put(i, idMap.get(mapBack.get(i)));
		g.idMap = newIdMap;
		return g;
	}

	public RtcGraph getOneNeighborhood (Integer c) {
		RtcGraph g = new RtcGraph();
		// We map the original integer IDs to new ones in order to minimize the range of vertex IDs
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapBack = new HashMap<Integer, Integer>();
		g.addVertex(0);										// add center vertex
		map.put(c, 0);
		mapBack.put(0, c);
		int index = 1;
		for (int v : adjacency.get(c)) {
			g.addVertex(index);
			map.put(v, index);
			mapBack.put(index, v);
			index++;
		}

		for (int v : g.nodes)
			for (int w : adjacency.get(mapBack.get(v)))
				if (map.containsKey(w) && map.get(w) > v)
					g.addEdge(v, map.get(w));

		// Construct idMap for new graph, mapping new integer IDs to original string identifiers
		HashMap<Integer, String> newIdMap = new HashMap<Integer, String>();
		for (int i = 0; i < g.size(); i++)
			newIdMap.put(i, idMap.get(mapBack.get(i)));
		g.idMap = newIdMap;
		return g;
	}

	public RtcGraph turnIntoSingleKernel () {
		RtcGraph g = new RtcGraph();
		// We map the original integer IDs to new ones in order to minimize the range of vertex IDs
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapBack = new HashMap<Integer, Integer>();

		// Copy vertices
		int index = 0;
		for (int v : nodes) {
			g.addVertex(index);
			map.put(v, index);
			mapBack.put(index, v);

			index++;
		}

		// Copy edges
		for (int v : g.nodes)
			for (int w : adjacency.get(mapBack.get(v)))
				if (map.get(w) > v)
					g.addEdge(v, map.get(w));

		// Construct idMap for new graph, mapping new integer IDs to original string identifiers
		HashMap<Integer, String> newIdMap = new HashMap<Integer, String>();
		for (int i = 0; i < g.size(); i++)
			newIdMap.put(i, idMap.get(mapBack.get(i)));
		g.idMap = newIdMap;
		return g;
	}

	public void deleteVertex (int vertex) {
		for (int neighbor : adjacency.get(vertex))
			adjacency.get(neighbor).remove(vertex);
		nodes.remove(vertex);
	}

	public void retainVertices (HashSet<String> vertices) {
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int v : nodes)
			if (!vertices.contains(getVertexName(v)))
				deleteSet.add(v);
		for (int v : deleteSet)
			deleteVertex(v);
	}

	public void retainVerticesByID (HashSet<Integer> vertices) {
		HashSet<Integer> deleteSet = new HashSet<Integer>();
		for (int v : nodes)
			if (!vertices.contains(v))
				deleteSet.add(v);
		for (int v : deleteSet)
			deleteVertex(v);
	}

	/**
	 * IMPORTANT:	UNDELETES ARE ONLY CORRECT when performed IN REVERSE ORDER of deletion!
	 *
	 * @param vertex
	 * @throws GraphException
	 */
	public void undeleteVertex (int vertex) {
		nodes.add(vertex);
		for (int neighbor : adjacency.get(vertex))
			adjacency.get(neighbor).add(vertex);
	}

	public RtcGraph getClone() {
		RtcGraph clone = new RtcGraph();
		for (int v : nodes)
			clone.addVertex(v);
		for (int v : nodes)
			for (int w : adjacency.get(v))
				if (w > v)
					clone.addEdge(v, w);
		clone.idMap = idMap;
		return clone;
	}

	public String getVertexName (Integer vertex) {
		return idMap.get(vertex);
	}

	public HashSet<String> getVertexNames (HashSet<Integer> ids) {
		HashSet<String> vertexNames = new HashSet<String>();
		for (int vertex : ids)
			vertexNames.add(idMap.get(vertex));
		return vertexNames;
	}

	public HashSet<String> getVertexNames () {
		HashSet<String> vertexNames = new HashSet<String>();
		for (int vertex : nodes)
			vertexNames.add(idMap.get(vertex));
		return vertexNames;
	}

	/** Returns the largest t, for which graph is a 2,t-club. t = n, if graph is a clique. 
	 * @throws GraphException */
	public int is2TClub (Model usedModel) {
		int t = this.size();
		for (int v : nodes)
		{
			for (int w : nodes)
			{
				if (w > v)
				{
					if (adjacent(v, w) && usedModel == Model.HEREDITARY)
						continue;

					int lengthAtMost2Paths = countCommonNeighbors(v,w);
					if (adjacent(v, w))
						lengthAtMost2Paths++;

					if (lengthAtMost2Paths == 0)
						return -1; // not a 2-club at all -> something is wrong

					if ((usedModel == Model.HEREDITARY || usedModel == Model.VB_MODEL) && lengthAtMost2Paths < t)
						t = lengthAtMost2Paths + 1;

					if (usedModel == Model.BICONNECTED)
					{
						int tmpT = noOfInternallyVertexDisjointPaths(v, w);
						if (tmpT < t)
							t = tmpT;
					}
				}
			}
		}
		return t;
	}

	/** Returns whether a given subset of vertices forms a (not necessarily maximal) clique. */
	public boolean isClique (HashSet<Integer> vertices) {
		for (int v : vertices) for (int w : vertices) if (v > w)
			if (!adjacent(v,w)) return false;
		return true;
	}

	/** Returns whether the whole graph is a clique. */
	public boolean isClique () {
		if (nodes.size() * (nodes.size() - 1) / 2 == getEdgeCount()) return true;
		else return false;
	}

	/** Returns whether the whole graph is connected. */
	public boolean isConnected ()
	{
		if (nodes.size() == 0)
			return true;
		int v = nodes.iterator().next();
		HashSet<Integer> visitedNodes = new HashSet<Integer>();
		Stack<Integer> activeNodes = new Stack<Integer>();
		activeNodes.add(v);
		visitedNodes.add(v);
		while (!activeNodes.isEmpty())
		{
			int currentNode = activeNodes.pop();
			for (int w : adjacency.get(currentNode))
			{
				if (!visitedNodes.contains(w))
				{
					activeNodes.add(w);
					visitedNodes.add(w);
				}
			}
		}
		return (visitedNodes.size() == nodes.size());
	}


	/**
	 * returns the edges of a path from v to w
	 * @param v
	 * @param w
	 * @return returns null if v=w or if there is no path. Otherwise the vertices of a path are returned.
	 */
	public LinkedList<Integer> getPath(int v, int w)
	{
		if (v == w)
		{
			return null;
		}
		long time = System.nanoTime();

		LinkedList<Integer> result = new LinkedList<Integer>();
		if (adjacent(v, w))
		{
			result.add(v);
			result.add(w);
			flowComputePath += System.nanoTime() - time;
			return result;
		}

		LinkedList<Integer> toVisit = new LinkedList<Integer>();
		HashMap<Integer,Integer> predecessor = new HashMap<Integer,Integer>();
		HashSet<Integer> seen = new HashSet<Integer>();

		toVisit.add(v);
		seen.add(v);

		boolean reached = false;

		while (!reached && !toVisit.isEmpty())
		{
			int activeVertex = toVisit.removeFirst();
			for (int neighbor : getNeighbors(activeVertex))
			{
				if (neighbor == w)
				{
					reached = true;
				}
				if (!seen.contains(neighbor))
				{
					predecessor.put(neighbor, activeVertex);
					seen.add(neighbor);
					toVisit.add(neighbor);
				}
			}
		}

		if (reached)
		{
			int activeVertex = w;
			while (activeVertex != v)
			{
				result.addFirst(activeVertex);
				activeVertex = predecessor.get(activeVertex);
			}
			result.addFirst(v);

			flowComputePath += System.nanoTime() - time;
			return result;
		}
		else
		{
			flowComputePath += System.nanoTime() - time;
			return null;
		}
	}

	// only for internal usage in maxflow computation where we use digraphs...
	private void augmentPath(LinkedList<Integer> path)
	{
		int lastVertex = path.removeFirst();
		while (!path.isEmpty())
		{
			int nextVertex = path.removeFirst();

			// reverse arc:
			adjacency.get(lastVertex).remove(nextVertex);
			adjacency.get(nextVertex).add(lastVertex);
			lastVertex = nextVertex;
		}
	}

	/**
	 * computes the number of (internally) vertex disjoint paths from v to w
	 * @param v
	 * @param w
	 * @return
	 */
	public int noOfInternallyVertexDisjointPaths(int v, int w)
	{
		return noOfInternallyVertexDisjointPaths(v,w,Integer.MAX_VALUE);
	}

	/**
	 * Checks whether at least {@code bound} the number of (internally) vertex disjoint paths from v to w
	 * @param v
	 * @param w
	 * @param bound
	 * @return
	 */
	public boolean enoughInternallyVertexDisjointPaths(int v, int w, int bound)
	{
		return noOfInternallyVertexDisjointPaths(v,w,bound) >= bound;
	}


	/**
	 * computes the number of (internally) vertex disjoint paths from v to w 
	 * @param v
	 * @param w
	 * @param upperBoundToStop 	When at least this number of (internally) vertex disjoint paths from v to w, then the function stops searching for more. 
	 * 							This improves the performance a bit. 
	 * @return
	 */
	private int noOfInternallyVertexDisjointPaths(int v, int w, int upperBoundToStop)
	{
		if (!nodes.contains(v) || !nodes.contains(w) || v == w)
			return 0;

		long time = System.nanoTime();

		// use standard reduction to the problem of finding edge-disjoint paths:
		// replace each vertex v by two vertices v1 and v2: one vertex v1 for the ingoing arcs and one vertex v2 for the outgoing arcs (plus one arc from v1 to v2). For simplicity we set v1 := v.
		int maxVertexNo = 0;
		for (int vertex : nodes)
		{
			if (maxVertexNo < vertex)
				maxVertexNo = vertex;
		}

		RtcGraph residualGraph = new RtcGraph();

		for (int vertex : nodes)
		{
			residualGraph.addVertex(vertex);
			int newVertex = vertex + maxVertexNo + 1;
			residualGraph.addVertex(newVertex);
			residualGraph.addEdge(vertex, newVertex);

			HashSet<Integer> neighbors = (HashSet<Integer>) this.getNeighbors(vertex).clone();

			HashSet<Integer> newNeighbors = new HashSet<Integer>();
			newNeighbors.add(newVertex);

			residualGraph.adjacency.put(vertex, newNeighbors); // we build a directed graph so we cannot use addEdge
			residualGraph.adjacency.put(newVertex, neighbors); // we build a directed graph so we cannot use addEdge
		}

		v = v + maxVertexNo + 1; // new starting vertex since v should be used in many paths

		flowConstructGraph += System.nanoTime() - time;

		int result = 0;
		LinkedList<Integer> path;
		while ( (path = residualGraph.getPath(v, w)) != null )
		{
			result++;

			if (upperBoundToStop <= result)
				return upperBoundToStop;

			residualGraph.augmentPath(path);
		}

		return result;
	}


	public void writeToFile(String directory, String filename, String[] header) {
		writeToFile(directory, filename, header, 0);
	}


	/** Writes this graph to disk in edge list format.
	 *
	 * @param directory String representation of the directory to which the graph shall be saved.
	 * @param filename	Name of the file to which the graph shall be saved.
	 * @param header	String array of header information for file. Pass null, if not required.
	 */
	public void writeToFile(String directory, String filename, String[] header, int indexDecrement) {
		File file; File dir;
		if (directory != null) {
			file = new File(directory+"/"+filename+".dimacs");
			dir = new File(directory+"/");
			if (!dir.exists()) dir.mkdir();
		} else {
			file = new File(filename+".dimacs");
		}
		try {
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			// Write header information to file
			if (header != null)
				for (int line = 0; line < header.length; line++)
					bw.write(header[line]+"\n");
			// Write adjacency information to file as edge list
			HashSet<Integer> seenNodes = new HashSet<Integer>();
			for (Integer v : nodes) {
				seenNodes.add(v);
				for (Integer w : adjacency.get(v))
					if (!seenNodes.contains(w))
						if (indexDecrement == 0)
							bw.write(idMap.get(v) + "\t" + idMap.get(w) + "\n");
						else {
							bw.write((Integer.parseInt(idMap.get(v)) - indexDecrement) + "\t" +
									(Integer.parseInt(idMap.get(w)) - indexDecrement) + "\n");
						}
			}
			bw.flush();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	public static RtcGraph ErdosRenyi (int n, double p) {
		RtcGraph g = new RtcGraph();
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		for (int i = 0; i < n; i++) {
			g.addVertex(i);
			map.put(i, Integer.toString(i));
		}
		for (int v = 0; v < n; v++)
			for (int w = v+1; w < n; w++)
				if (Math.random() < p)
					g.addEdge(v, w);
		g.idMap = map;
		return g;
	}

	public static RtcGraph Gendreau (int n, double a, double b) {
		RtcGraph g = new RtcGraph();
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		double[] probabilities = new double[n];
		for (int i = 0; i < n; i++) {
			g.addVertex(i);
			map.put(i, Integer.toString(i));
			probabilities[i] = Math.random()*(b-a) + a;
		}
		for (int v = 0; v < n; v++)
			for (int w = v+1; w < n; w++)
				if (Math.random() < (probabilities[v] + probabilities[w])/2)
					g.addEdge(v, w);
		g.idMap = map;
		return g;
	}

	public static RtcGraph AlternateGendreau (int n, double a, double b) {
		RtcGraph g = new RtcGraph();
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		double[] probabilities = new double[n];
		for (int i = 0; i < n; i++) {
			g.addVertex(i);
			map.put(i, Integer.toString(i));
			probabilities[i] = Math.random()*(b-a) + a;
		}
		for (int v = 0; v < n; v++)
			for (int w = v+1; w < n; w++)
				if (Math.random() < probabilities[v] * probabilities[w])
					g.addEdge(v, w);
		g.idMap = map;
		return g;
	}

	/**
	 *
	 * @param n
	 * @param deg
	 * @param seeds
	 * @param cycle		Set to true to seed from a cycle and false to seed from a clique.
	 * @return
	 * @throws GraphException
	 */
	public static RtcGraph BarabasiAlbert (int n, int deg, int seeds, boolean cycle) {
		RtcGraph g = new RtcGraph();
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		int degreeSum = (seeds - 1) * seeds;
		if (cycle) degreeSum = 2 * seeds;
		int[] accumulatedDegrees = new int[n];
		if (cycle) for (int i = 0; i < seeds; i++)
			accumulatedDegrees[i] = (i+1) * 2;
		if (!cycle) for (int i = 0; i < seeds; i++)
			accumulatedDegrees[i] = (i+1) * (seeds - 1);
		for (int i = 0; i < seeds; i++) {
			g.addVertex(i);
			map.put(i, Integer.toString(i));
			if (cycle) {
				if (i>0) g.addEdge(i, i-1);
				if (i==seeds-1) g.addEdge(i, 0);
			} else for (int j = 0; j < i; j++)
				g.addEdge(i, j);
		}
		for (int i = seeds; i < n; i++) {
			g.addVertex(i);
			map.put(i, Integer.toString(i));
			for (int e = 0; e < deg; e++) {
				int pickedVertex = -1;
				while (pickedVertex == -1 || g.adjacent(i, pickedVertex)) {
					int rand = (int) (Math.random() * degreeSum);
					int vertexPointer = 0;
					while (accumulatedDegrees[vertexPointer] < rand)
						vertexPointer++;
					pickedVertex = vertexPointer;
				}
				g.addEdge(i, pickedVertex);
				degreeSum += 1;
				for (int j = pickedVertex; j < i; j++)
					accumulatedDegrees[j]++;
			}
			degreeSum += deg;
			accumulatedDegrees[i] = accumulatedDegrees[i-1] + deg;
		}
		g.idMap = map;
		return g;
	}

	public static RtcGraph RandomRoughlyRegular (int n, int deg) {
		if (n < deg-1 || n % 2 + deg % 2 == 2) return null;
		RtcGraph g = new RtcGraph();
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		HashMap<Integer, ArrayList<Integer>> degreeBuckets = new HashMap<Integer, ArrayList<Integer>>();
		for (int d = 0; d <= deg+1; d++)
			degreeBuckets.put(d, new ArrayList<Integer>());
		for (int v = 0; v < n; v++) {
			g.addVertex(v);
			map.put(v, Integer.toString(v));
			degreeBuckets.get(0).add(v);
		}
		int edges = 0;
		while (edges != n*deg/2) {
			int bucketPtr = 0;
			while (degreeBuckets.get(bucketPtr).size() == 0)
				bucketPtr++;
			int firstVertexPos = (int)(Math.random() * (degreeBuckets.get(bucketPtr).size()));
			int firstVertex = degreeBuckets.get(bucketPtr).remove(firstVertexPos);
			int bucketPtr2 = bucketPtr;
			while (degreeBuckets.get(bucketPtr2).size() == 0)
				bucketPtr2++;
			while (bucketPtr2 <= deg) {
				boolean stop = false;
				for (int v : degreeBuckets.get(bucketPtr2)) {
					if (!g.adjacent(firstVertex, v)) {
						stop = true; break;
					}
				}
				if (stop) break;
				bucketPtr2++;
			}
			int secondVertexPos = (int)(Math.random() * (degreeBuckets.get(bucketPtr2).size()));
			int secondVertex = degreeBuckets.get(bucketPtr2).remove(secondVertexPos);
			if (!g.adjacent(firstVertex, secondVertex)) {
				g.addEdge(firstVertex, secondVertex);
				degreeBuckets.get(bucketPtr+1).add(firstVertex);
				degreeBuckets.get(bucketPtr2+1).add(secondVertex);
				edges++;
			} else {
				degreeBuckets.get(bucketPtr).add(firstVertex);
				degreeBuckets.get(bucketPtr2).add(secondVertex);
			}
		}
		g.idMap = map;
		return g;
	}

	public RtcGraph () {
		idMap = new HashMap<Integer, String>();
		nodes = new HashSet<Integer>();
		adjacency = new HashMap<Integer, HashSet<Integer>>();
	}

	public void addVertex (Integer v) {
		nodes.add(v);
		adjacency.put(v, new HashSet<Integer>());
	}

	public void addVertex (Integer v, String name) {
		nodes.add(v);
		adjacency.put(v, new HashSet<Integer>());
		idMap.put(v, name);
	}


	public void addEdge (Integer v, Integer w) {
		if (v == w) return;		// No loops!
		adjacency.get(v).add(w);
		adjacency.get(w).add(v);
	}

	public void deleteEdge (Integer v, Integer w) {
		if (v == w) return;		// No loops!
		adjacency.get(v).remove(w);
		adjacency.get(w).remove(v);
	}

	public HashMap<String, Integer> getReversedIdMap(){
		HashMap<String, Integer> map = new HashMap<>();

		for(Integer v : nodes)
			map.put(idMap.get(v), v);

		return map;
	}

	public int getHighestVertexID(){
		int max = 0;
		for(int id : nodes)
			if(id > max)
				max = id;
		return max;
	}

	private void buildGraphFromFile(String filename, int type){
		HashMap<String, Integer> seenIds = new HashMap<String, Integer>();	// only required during construction
		File file = new File(filename);

		try {
			final LineNumberReader reader = new LineNumberReader(new FileReader(file));
			String str;
			int id = 0;
			if (type == 0) {		// Parse DIMACS file
				reader.readLine();  //read first line containing the number of vertices and edges
				while ((str = reader.readLine()) != null)
				{
					str = str.trim();										// trim away whitespace at either end of line
					if (!str.startsWith("#") && !str.startsWith("%")) {								// skip comment lines
						StringTokenizer tokens = new StringTokenizer(str);
						if (tokens != null && tokens.countTokens() > 1) {	// only consider well-formed lines
							String vertexA = tokens.nextToken();
							String vertexB = tokens.nextToken();
							if (!seenIds.containsKey(vertexA)) {		// add vertex 0 if never seen before
								seenIds.put(vertexA, id);
								addVertex(id, vertexA);
								id++;
							}
							if (!seenIds.containsKey(vertexB)) {		// add vertex 1 if never seen before
								seenIds.put(vertexB, id);
								addVertex(id, vertexB);
								id++;
							}
							// add edge to both adjacency lists
							addEdge(seenIds.get(vertexA), seenIds.get(vertexB));
						}
					}
				}
			} else if (type == 1) {		// parse METIS file
				StringTokenizer header = new StringTokenizer(reader.readLine());
				int vertices = Integer.parseInt(header.nextToken());
				int edges = Integer.parseInt(header.nextToken());
				boolean weights = header.hasMoreTokens() && !header.nextToken().equals("0");
				for (int v = 0; v < vertices; v++) {
					addVertex(v, Integer.toString(v+1));
				}
				int edgesCounted = 0;
				int currentVertex = 0;
				while ((str = reader.readLine()) != null)
				{
					if (!str.startsWith("#") && !str.startsWith("%")) {				// skip comment lines
						StringTokenizer tokens = new StringTokenizer(str);
						while (tokens.hasMoreTokens()) {
							int v = Integer.parseInt(tokens.nextToken()) - 1;
							if (!adjacent(currentVertex, v)) {
								addEdge(currentVertex, v);
								edgesCounted++;
							}
							if (weights) tokens.nextToken(); // Throw away weight token
						}
						currentVertex++;
						tokens = null;
					}
				}
				if (edges != edgesCounted) {
					System.out.println("WARNING: Wrong edge count while loading METIS file "+filename+".");
					System.out.println("Edges stated: "+edges+"   Edges counted: "+edgesCounted+"\n");
				}
			} else if (type == 2) {
				while ((str = reader.readLine()) != null)
				{
					StringTokenizer tokens = new StringTokenizer(str);
					String lineType = tokens.nextToken();
					if (lineType.equals("c")) continue;
					if (lineType.equals("p")) {
						tokens.nextToken();
						int vertices = Integer.parseInt(tokens.nextToken());
						for (int v = 0; v < vertices; v++) {
							addVertex(v, Integer.toString(v+1));
						}
					}
					if (lineType.equals("e")) {
						addEdge(Integer.parseInt(tokens.nextToken()) - 1, Integer.parseInt(tokens.nextToken()) - 1);
					}
				}
			} else if (type == 3) {		// Edges only
				while ((str = reader.readLine()) != null)
				{
					str = str.trim();										// trim away whitespace at either end of line
					if (!str.startsWith("#") && !str.startsWith("%")) {								// skip comment lines
						StringTokenizer tokens = new StringTokenizer(str);
						if (tokens != null && tokens.countTokens() > 1) {	// only consider well-formed lines
							String vertexA = tokens.nextToken();
							String vertexB = tokens.nextToken();
							if (!seenIds.containsKey(vertexA)) {		// add vertex 0 if never seen before
								seenIds.put(vertexA, id);
								addVertex(id, vertexA);
								id++;
							}
							if (!seenIds.containsKey(vertexB)) {		// add vertex 1 if never seen before
								seenIds.put(vertexB, id);
								addVertex(id, vertexB);
								id++;
							}
							// add edge to both adjacency lists
							addEdge(seenIds.get(vertexA), seenIds.get(vertexB));
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			System.out.println("Could not locate input file '"+filename+"'.");
			System.exit(0);
		}
	}

	/**
	 * Creates a subgraph of this graph that only contains the specified vertices/nodes. The given set of vertices is
	 * used as the node set (not a copy).
	 * @param vertices a subset of nodes
	 * @return a new graph with only the specified vertices
	 */
	public RtcGraph getSubgraph(HashSet<Integer> vertices){
		RtcGraph solution = new RtcGraph();
		solution.nodes = vertices;
		solution.idMap = this.idMap;
		solution.adjacency = this.getLocalAdjacency(vertices);
		return solution;
	}

	/**
	 * Computes N_2(v) \cap N_2(u) for two vertices v, u
	 * @param v first vertex id
	 * @param u second vertex id
	 * @return a new HashSet with vertices of to N_2(v) \cap N_2(u)
	 */
	public HashSet<Integer> getCommonTwoNeighborhood(int v, int u){
		HashSet<Integer> N2v = new HashSet<>(), N2u = new HashSet<>();

		getNeighbors(v).forEach(x -> {
			N2v.addAll(getNeighbors(x));
			N2v.add(x);
		});

		getNeighbors(u).forEach(x -> {
			if(N2v.contains(x))
				N2u.add(x);
			for(Integer y : getNeighbors(x))
				if(N2v.contains(y))
					N2u.add(y);
		});

		return N2u;
	}

	/**
	 * Creates a copy of the adjacency for a specified neighborhood.
	 * @param neighborhood a subset of vertices of this graph
	 * @return a deep copy of the local adjacency
	 */
	public HashMap<Integer, HashSet<Integer>> getLocalAdjacency(HashSet<Integer> neighborhood){
		HashMap<Integer, HashSet<Integer>> local = new HashMap<>();
		HashSet<Integer> neighbors;
		for(int v : neighborhood){
			if(!nodes.contains(v))
				continue;
			neighbors = new HashSet<>();
			for(int w : getNeighbors(v))
				if(neighborhood.contains(w))
					neighbors.add(w);
			local.put(v, neighbors);
		}
		return local;
	}

	/**
	 * Computes the density of the graph.
	 * @return the density
	 */
	public float getDensity(){
		int edges = 0;
		for(int v : nodes){
			edges += adjacency.get(v).size();
		}
		return (float) edges / (float) (nodes.size() * nodes.size() - 1);
	}
}
