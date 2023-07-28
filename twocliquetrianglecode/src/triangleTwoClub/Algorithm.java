package triangleTwoClub;

import robustTwoClub.graph.Edge;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.dataReduction.*;
import triangleTwoClub.graph.T2CGraph;

import java.util.*;

public class Algorithm implements ITriangleTwoCliqueAlgorithm {

    private static boolean checkSolution = true;
    private static float TNR_LCR_densityThreashold = 0.05f;

    private ConditionChecker conditionChecker;
    private Output output;
    private int ell, highestVertexID;
    private int[] startPointSizes;
    private T2CGraph graph;

    private T2CGraph kopie;

    private RtcGraph original;
    private RtcGraph inputGraph, currentBestSolution;
    private HashMap<Integer, HashSet<Triangle>> triangleMap;
    private boolean useEdgeCondition, abortCurrentBranch, abortCurrentRun, timeLimitEnabled, timeLimitReached,
            printInformationToConsole, branchingAborted, addSolutionToOutput;
    private ArrayList<DataReductionRule> dataReductionRules;
    private long currentRunStartTime, algorithmTimeLimit, abortTime;
    private UpperBound upperBound;
    private LowerBound.LowerBoundVariant lowerBoundVariant;
    private T2CGraph.ConflictGraphType conflictGraphType;

    public Algorithm(){
        ell = 1;
        useEdgeCondition = false;
        timeLimitEnabled = false;
        addSolutionToOutput = false;
        algorithmTimeLimit = 0L;
        lowerBoundVariant = LowerBound.LowerBoundVariant.LB1;
        conflictGraphType = T2CGraph.ConflictGraphType.SIMPLE;
        resetAlgorithm();
        initDataReductionRules();
    }

    public RtcGraph run(RtcGraph inputGraph, String graphName, String outputFile, boolean printInformation){
        // prepare
        printInformationToConsole = printInformation;
        print("starting ...");
        currentBestSolution = null;
        this.inputGraph = inputGraph;
        prepareOutput(graphName);
        this.highestVertexID = inputGraph.getHighestVertexID();
        HashSet<Integer> initialVertices = (HashSet<Integer>)inputGraph.getVertices().clone();
        graph = new T2CGraph(inputGraph, this);
        kopie = new T2CGraph(inputGraph, this);
        graph.setConflictType(conflictGraphType);
        kopie.setConflictType(conflictGraphType);
        LowerBound lowerBound = new LowerBound(this, lowerBoundVariant);
        upperBound = new UpperBound(this);

        this.original = inputGraph.getClone();

        // run
        currentRunStartTime = System.currentTimeMillis();
        abortTime = currentRunStartTime + algorithmTimeLimit;
        new VertexDegreeRule(this, "InitialDeletionRule").apply();
        computeTriangleMap();
        if(triangleMap != null){
            // triangleMap is null if the triangle computation was aborted by the time limit
            doInitialDataReduction();
            removeNonExistingTriangles(initialVertices);
        }
        print("searching for lower bound solution ...");
        long lowerBoundTimeLimit = -1;
        if(timeLimitEnabled)
            lowerBoundTimeLimit = abortTime <= System.currentTimeMillis() ? 0 : abortTime - System.currentTimeMillis();
        lowerBound.setTimeLimit(lowerBoundTimeLimit);
        currentBestSolution = lowerBound.getLowerBoundSolution(); // compute lower bound
        print("lower bound = " + getBestSolutionSize());
        if(lowerBoundVariant == LowerBound.LowerBoundVariant.LB1
                || (lowerBoundVariant == LowerBound.LowerBoundVariant.LB2 && currentBestSolution != null
                && currentBestSolution.size() > output.getLoggedValue(LowerBound.LB1SizeForLB2)))
            new TwoNeighborhoodRule(this,
                    currentBestSolution == null ? 0 : currentBestSolution.size()).apply();
        applyDeletionRules();
        chooseTnrOrLcrRule();
        ArrayList<Integer> startPoints = prepareStartPoints();
        output.log("StartPointCount", startPoints.size());
        print("number of start points = " + startPoints.size());

        print("start branching ...");
        // Track time needed for preprocessing
        output.log("PreProcTime[s]", System.currentTimeMillis() - currentRunStartTime);
        for(int i = 0; i < startPoints.size(); i++){
            if(abortCurrentRun) break;
            output.log("CurrentStartPoint", 1);
            searchAtStartPoint(startPoints.get(i), kopie);
        }

        output.log("Time[s]", System.currentTimeMillis() - currentRunStartTime); // log algorithm runtime


        // clean up and validate the result
        int n = -1, m = -1;
        if(currentBestSolution != null){
            n = currentBestSolution.size();
            m = currentBestSolution.getEdgeCount();
        }
        applyDataReductionToGraph(currentBestSolution, useEdgeCondition, ell);
        if(currentBestSolution != null && (n != currentBestSolution.size() || m != currentBestSolution.getEdgeCount())){
            System.out.print("Warning: Something was removed from the solution! ");
            System.out.print("(Graph:  " + graphName + ", L=" + ell + ", edgeCondition=" + useEdgeCondition + ") ");
            System.out.println("n: " + n + " -> " + currentBestSolution.size() + " , m: " + m + " -> "
                    + currentBestSolution.getEdgeCount());
        }
        if(useEdgeCondition && currentBestSolution != null && checkSolution
                && !ConditionChecker.isEdgeLTriangleTwoClub(currentBestSolution,
                Triangle.getTriangleHashMap(currentBestSolution), ell)){
            System.out.print("Warning: Solution is not an edge triangle 2 club! ");
            if(!printInformation)
                System.out.print("Configuration: " + graphName + ", L=" + ell + ", edge=" + useEdgeCondition + ")");
            System.out.println();
        }
        if(!useEdgeCondition && currentBestSolution != null && checkSolution
                && !ConditionChecker.isLTriangleTwoClub(currentBestSolution,
                Triangle.getTriangleHashMap(currentBestSolution), ell)){
            System.out.print("Warning: Solution is not an vertex triangle 2 club!");
            if(!printInformation)
                System.out.print(" (Configuration: " + graphName + ", L=" + ell + ", edge=" + useEdgeCondition + ")");
            System.out.println();
        }

        // finish
        output.addInformation("BranchingAborted", String.valueOf(branchingAborted));
        output.log("SolutionSize", currentBestSolution == null ? 0 : currentBestSolution.size());
        boolean timeLimit = timeLimitEnabled && (timeLimitReached
                || output.getInformation(LowerBound.lowerBoundTimeLimitAbortName).equals("true")
                || output.getInformation("TriangleComputationAborted").equals("true"));
        output.addInformation("TimeLimitReached", "" + timeLimit);
        if(addSolutionToOutput)
            output.logGraph("Solution", currentBestSolution);
        if(outputFile != null)
            output.writeOutputFile(outputFile);
        if(timeLimitReached) print("algorithm stopped, time limit reached!");
        print("finished after " + output.getLoggedValue("Time[s]").intValue() + "ms");
        print("result size = " + getBestSolutionSize());
        resetAlgorithm();
        return currentBestSolution;
    }

    /**
     * Searches for a l-triangle-2-club at the given start point
     * @param startPointID the start vertex id
     */
    private void searchAtStartPoint(int startPointID, T2CGraph kopie){
        if(abortCurrentRun = (timeLimitEnabled && System.currentTimeMillis() > abortTime)){
            timeLimitReached = true;
            branchingAborted = true;
            return;
        }

        T2CGraph globalGraph = kopie; // save global graph for later // input graph
        if(isStartPointValid(startPointID)){
            HashSet<Integer> localVertices = graph.getGraph().getTwoNeighbors(startPointID);
            localVertices.add(startPointID);
            graph = new T2CGraph(globalGraph.getGraph().getSubgraph(localVertices),this);
            graph.markVertex(startPointID);
            applyDeletionRules();
            graph.buildConflictGraph();
            abortCurrentBranch = false;

            /**Map<Integer, Set<Integer>> 	neighbours =  identifyVerticesWithDistanceMaxTwo(graph);
            graph.deleteVertex(6);
            Set<Integer> vertex = graph.getGraph().getVertices();
            Set<Set<Integer>> setOfTriangle = getSetsOfTriangle();
            Set<Integer> clique = getClique(startPointID,setOfTriangle,neighbours);


            if(getBestSolutionSize() < clique.size()){ // wenn Clique größer doese nehmen,
                Set<Integer> eliminationList = getElementsOfSetOneNotInSetTwo(graph.getGraph().getVertices(),clique);
                RtcGraph temp = graph.getGraph().getClone();
                for (Integer toDeleat : eliminationList) {
                    temp.deleteVertex(toDeleat);
                }
                currentBestSolution = temp;
            }*/
            branchClique(1);
            //branch(1);
        }
        graph = globalGraph;
        graph.deleteVertex(startPointID);
    }

    private void branch(int depth){
        // algorithm management
        output.log("BranchCount", 1);
        if(abortCurrentRun) return;
        if(output.getLoggedValue("MaxBranchDepth") < depth)
            output.log("MaxBranchDepth", depth, false);
        if(printInformationToConsole)
            output.printStatusToConsole();
        if(timeLimitEnabled && System.currentTimeMillis() > abortTime){
            abortCurrentRun = true;
            timeLimitReached = true;
            branchingAborted = true;
            return;
        }

        // apply data reduction
        Edge rollbackPoint = graph.getRollbackPoint();
        applyDataReductionRules();
        if(abortCurrentBranch) {
            abortCurrentBranch = false;
            graph.rollback(rollbackPoint);
            return;
        }

        // check current branching state is a solution
        if(useEdgeCondition ? conditionChecker.isEdgeLTriangleTwoClub() : conditionChecker.isLTriangleTwoClub()){ // darf auf keinen Fall auskommentiert werden
            if(getBestSolutionSize() < graph.getGraph().size()){
                // this solution is the new best solution
                currentBestSolution = graph.getGraph().getClone();
            }
            return;
        }

        // check if the graph is big enough to continue branching
        if(graph.getGraph().size() <= getBestSolutionSize())
            return;

        // find a vertex to continue branching
        int v = conditionChecker.getMostIncompatibleVertexID(inputGraph); //hier wäre es gut den Inputgraph zu übergeben, aber das klappt nicht
        Edge branchingRollbackPoint;



        // option 1: branch without v
        branchingRollbackPoint = graph.getRollbackPoint();
        graph.deleteVertex(v);
        branch(depth + 1);
        graph.rollback(branchingRollbackPoint);

        // option 2: branch with v marked
        branchingRollbackPoint = graph.getRollbackPoint();
        graph.markVertex(v);
        branch(depth + 1);
        graph.rollback(branchingRollbackPoint);

        graph.rollback(rollbackPoint); // rollback to the state before data reduction
    }

    private void branchClique(int depth){
        // algorithm management
        output.log("BranchCount", 1);
        if(abortCurrentRun) return;
        if(output.getLoggedValue("MaxBranchDepth") < depth)
            output.log("MaxBranchDepth", depth, false);
        if(printInformationToConsole)
            output.printStatusToConsole();
        if(timeLimitEnabled && System.currentTimeMillis() > abortTime){
            abortCurrentRun = true;
            timeLimitReached = true;
            branchingAborted = true;
            return;
        }

        // apply data reduction
        Edge rollbackPoint = graph.getRollbackPoint();
        applyDataReductionRules();
        if(abortCurrentBranch) {
            abortCurrentBranch = false;
            graph.rollback(rollbackPoint);
            return;
        }

        // check current branching state is a solution
        if( conditionChecker.isLTriangleTwoClique()){ // darf auf keinen Fall auskommentiert werden
            if(getBestSolutionSize() < graph.getGraph().size()){
                // this solution is the new best solution
                currentBestSolution = graph.getGraph().getClone();
            }
            return;
        }

        // check if the graph is big enough to continue branching
        if(graph.getGraph().size() <= getBestSolutionSize())
            return;

        if(graph.getGraph().size() <= 2)
            return;

        int maxDepth = (int)Math.pow(getOriginalGraph().getGraph().size(),2);
        if(depth > maxDepth )
            return;

        // find a vertex to continue branching
        int v = conditionChecker.getMostIncompatibleVertexIDClique(); //hier wäre es gut den Inputgraph zu übergeben, aber das klappt nicht
        Edge branchingRollbackPoint;



        // option 1: branch without v
        branchingRollbackPoint = graph.getRollbackPoint();
        graph.deleteVertex(v);
        branchClique(depth + 1);
        graph.rollback(branchingRollbackPoint);

        // option 2: branch with v marked
        branchingRollbackPoint = graph.getRollbackPoint();
        //if(branchingRollbackPoint != null){
        graph.markVertex(v);
        branchClique(depth + 1);
        graph.rollback(branchingRollbackPoint);
        //}


        graph.rollback(rollbackPoint); // rollback to the state before data reduction
    }

    //Todo Beurteilung von Algorithmus zur erstellung einer Map von Kanten

    public Map<Integer, Set<Integer>> identifyVerticesWithDistanceMaxTwo(T2CGraph graph){
        RollbackRtcGraph internGraph = graph.getGraph();
        Set<Integer> vertices = internGraph.getVertices();
        Map<Integer,Set<Integer>> output = new HashMap<>();
        for (Integer v : vertices){
            output.put(v,internGraph.getTwoNeighbors(v));
        }
        return (output);
    }

    public Set<Set<Integer>> getSetsOfTriangle(){
        Set<Set<Integer>> set = new HashSet<>();
        for (int i : triangleMap.keySet()){
            for (Triangle triangle : triangleMap.get(i)){
                set.add(Set.copyOf(triangle.getVertices()));
            }
        }
        return (set);
    }


    @Override
    public Map<Integer, Set<Set<Integer>>> mapRemainingTrianglesInGraph(T2CGraph globalGraph){
        Set<Integer> vertex = globalGraph.getGraph().getVertices();
        Set<Set<Integer>> setTriangle = getSetsOfTriangle();
        Map<Integer, Set<Set<Integer>>> output = new HashMap<>();

        for (Set<Integer> set : setTriangle){
            if(overlapOfTwoSets(set,vertex).containsAll(set)){
                output = addTriangleToMap(output,set);
            }
        }

        return output;
    }

    @Override
    public boolean containsOnlyTriangle(RtcGraph graph) {

        Set<Set<Integer>> setsOfTriangle = getSetsOfTriangle();
        Set<Integer> vertices = graph.getVertices();
        Set<Integer> usedVertices = new HashSet<>();
        for (Set<Integer> triangle : setsOfTriangle){
            if(vertices.containsAll(triangle)){
                usedVertices.addAll(triangle);
            }
        }
        return (usedVertices.containsAll(vertices));

    }

    @Override
    public boolean isClique(RtcGraph graph) {
        Set<Integer> vertices = graph.getVertices();
        for (Integer v : vertices){
            Set<Integer> neighbors = new HashSet<>(getOriginalGraph().getGraph().getTwoNeighbors(v));
            neighbors.add(v);
            if (!neighbors.containsAll(vertices)){
                return(false);
            }
        }
        return (true);
    }

    private Map<Integer, Set<Set<Integer>>> addTriangleToMap(Map<Integer,Set<Set<Integer>>> map,Set<Integer> triangle){
        for (Integer i : triangle){
            if(map.containsKey(i)){
                map.get(i).add(triangle);
            }else{
                Set<Set<Integer>> temp = new HashSet<>();
                temp.add(triangle);
                map.put(i,temp);
            }
        }
        return (map);
    }

    public Set getElementsOfSetOneNotInSetTwo(Set setOne, Set setTwo){
        Set tempSetOne = new HashSet(setOne);

        tempSetOne.removeAll(setTwo);
        return (tempSetOne);
    }

    public Set overlapOfTwoSets(Set setA, Set setB){
        Set tempSet = new HashSet(setA);

        Set removalSetA = getElementsOfSetOneNotInSetTwo(setA,setB);
        Set removalSetB = getElementsOfSetOneNotInSetTwo(setB,setA);

        tempSet.removeAll(removalSetA);
        tempSet.removeAll(removalSetB);

        return tempSet;
    }

    @Override
    public T2CGraph getOriginalGraph() {
        return new T2CGraph(original, this);
    }


    public Set<Integer> getClique(Integer startPoint,Set<Set<Integer>> setsOfTriangle,Map<Integer, Set<Integer>> mapOfVerticesInRange){
        Set<Integer> output = new HashSet<>();
        if (!mapOfVerticesInRange.containsKey(startPoint)){
            return (output);
        }

        Set<Integer> setPossibleVertices = new HashSet<>();
        for (Set<Integer> triangle : setsOfTriangle){
            setPossibleVertices.addAll(triangle);
        }

        output = new HashSet<>(mapOfVerticesInRange.get(startPoint));
        Set<Integer> notInAnyTriangle = getElementsOfSetOneNotInSetTwo(output,setPossibleVertices);

        output.removeAll(notInAnyTriangle);
        Set<Integer> copyVertices = new HashSet<>(output);
        for (Integer vertices : copyVertices){
            Set<Integer> reachebleVertices = mapOfVerticesInRange.get(vertices);
            reachebleVertices.add(vertices);
            output = overlapOfTwoSets(output,reachebleVertices);
        }
        output.add(startPoint);

        return(output);

    }

    /**
     * Sets parameter l.
     * @param l the number of triangles each vertex/edge has to be part of
     */
    public void setL(int l){
        if(l > 0)
            this.ell = l;
        else
            System.out.println("L must be at least 1");
    }

    /**
     * Enables/disables the edge condition.
     * If it is enabled each edge has to be part of at leas l triangles. Otherwise only vertices have to be part
     * of at least l triangles.
     * @param b the value for the edge condition
     */
    public void enableEdgeCondition(boolean b){
        this.useEdgeCondition = b;
        initDataReductionRules(); // re-initialize with edge rule
    }

    /**
     * Sets an time limit to the algorithm. A value of zero or lower will disable the time limit.
     * When the time limit is reached, the algorithm stops at the next start of a branch.
     * @param timeLimitMilliSeconds the time limit in milliseconds
     */
    public void setTimeLimit(long timeLimitMilliSeconds){
        if(timeLimitMilliSeconds > 0){
            timeLimitEnabled = true;
            algorithmTimeLimit = timeLimitMilliSeconds;
        }else{
            timeLimitEnabled = false;
        }
    }

    /**
     * Sets the lower bound the algorithm uses.
     * @param variant a variant of the lower bound
     */
    public void setLowerBoundVariant(LowerBound.LowerBoundVariant variant){
        this.lowerBoundVariant = variant;
    }

    /**
     * Specifies how the conflict graph is built.
     * @param t the type of conflict graph
     */
    public void setConflictGraphType(T2CGraph.ConflictGraphType t){
        this.conflictGraphType = t;
        initDataReductionRules(); // re-init without maximum matching rule if DISABLED
    }

    /**
     * Enables or disables the output of the solution graph to the other output.
     * @param flag true if you want the solution in the output
     */
    public void addSolutionToOutput(boolean flag){
        this.addSolutionToOutput = flag;
    }

    /**
     * Checks if a start point can produce a better solution than the current best one.
     * @param startPointID the id of the start vertex
     * @return true if this start point may produce a better solution
     */
    private boolean isStartPointValid(int startPointID){ //TODO überprüfe ob diese Fubnktion wichtig und richtig ist
        // check if start point can produce a better solution
        if(currentBestSolution != null){
            if(upperBound.getUpperBound(startPointID) <= getBestSolutionSize()){
                output.log("UpperBoundSkips", 1);
                return false;
            }
            // check if neighborhood is large enough
            if(getGraph().getTwoNeighbors(startPointID).size() + 1 <= getBestSolutionSize()){
                output.log("StartPointSkips", 1);
                return false;
            }
        }
        return true;
    }

    /**
     * Initializes the data reduction rules. The order of the rules in the list specifies the order of execution.
     */
    private void initDataReductionRules(){
        dataReductionRules = new ArrayList<>();

        // rules for vertex deletion
        /**    dataReductionRules.add(new VertexDegreeRule(this));
         dataReductionRules.add(new IncompatibleResolutionRule(this));
         if(useEdgeCondition)
         dataReductionRules.add(new EdgeTriangleRule(this)); */
        dataReductionRules.add(new VertexTriangleRule(this));
        /**     DataReductionRule d = new LowCompatibilityRule(this);
         d.setActivated(false); // deactivate LCR before branching, chooseTnrOrLcrRule() will add a new activated one
         dataReductionRules.add(d);
         //dataReductionRules.add(new NeighborhoodSupersetRule(this));

         // rules for marking vertices
         dataReductionRules.add(new NoChoiceRule2(this));
         dataReductionRules.add(new NoChoiceRule3(this));

         // rules to abort branches
         dataReductionRules.add(new MarkedIncompatibleRule(this));
         if(conflictGraphType != T2CGraph.ConflictGraphType.DISABLED)
         dataReductionRules.add(new MaximumMatchingRule(this)); */ // this rule requires a conflict graph
    }

    /**
     * Applies all data reduction rules. If after any rule the branch can be aborted, this method returns without
     * applying the remaining rules.
     */
    private void applyDataReductionRules(){
        applyDeletionRules();
        for(DataReductionRule d : dataReductionRules){
            if(d instanceof DeleteRule) continue; // DeleteRules were already applied beforehand
            d.apply();
            if(this.abortCurrentBranch)
                break;
        }
    }

    /**
     * Applies all data reduction rules that are of type DeleteRule.
     */
    private void applyDeletionRules(){
        /**  boolean verticesDeleted;

         outerLoop:
         do{
         verticesDeleted = false;
         for(DataReductionRule d : dataReductionRules){
         // skip other rules
         if(!(d instanceof DeleteRule)) continue;

         Double oldCount = output.getLoggedValue(d.getRuleName());
         d.apply();

         // check if branch can get aborted
         if(abortCurrentBranch)
         break outerLoop;

         // check if vertices were deleted
         if(!output.getLoggedValue(d.getRuleName()).equals(oldCount))
         verticesDeleted = true;

         }
         }while(verticesDeleted); */
    }

    /**
     * Only applies data reduction rules that are meant to be used before branching.
     */
    private void doInitialDataReduction(){ //lässt sich nicht auskommentieren
        for(DataReductionRule d : dataReductionRules)
            if(d.applyInitially())
                d.apply();
    }

    /**
     * TNR Rule is better than LCR on graphs with a small density.
     * Choose a rule based on the current graphs density.
     */
    private void chooseTnrOrLcrRule(){
        // find rule in reduction rules
        int ruleIndex;
        DataReductionRule r;
        try {
            ruleIndex = useEdgeCondition ? 4 : 3;
            r = dataReductionRules.get(useEdgeCondition ? 4 : 3);
        }catch(IndexOutOfBoundsException e){
            r = null;
            ruleIndex = -1;
        }
        if(r == null || (!(r instanceof LowCompatibilityRule) && !(r instanceof TwoNeighborhoodRule))){
            DataReductionRule rule;
            for(int i = 0; i < dataReductionRules.size(); i++){
                rule = dataReductionRules.get(i);
                if(rule instanceof LowCompatibilityRule || rule instanceof TwoNeighborhoodRule){
                    ruleIndex = i;
                    break;
                }
            }
        }

        if(ruleIndex < 0)
            return; // nothing to do

        dataReductionRules.remove(ruleIndex);
        String usedRule;
        if(getGraph().getDensity() > TNR_LCR_densityThreashold) {
            dataReductionRules.add(ruleIndex, new LowCompatibilityRule(this));
            usedRule = "LCR";
        }else{
            dataReductionRules.add(ruleIndex, new TwoNeighborhoodRule(this,
                    currentBestSolution == null ? 0 : currentBestSolution.size()));
            usedRule = "TNR";
        }
        output.addInformation("LCR/TNR", usedRule);
    }

    /**
     * Removes the triangle of removed vertices, a removed vertex is part of the inputGraph but not of the workingGraph.
     * This aims to reduce the number of triangles in sparse graphs.
     */
    private void removeNonExistingTriangles(HashSet<Integer> initialVertices){
        for(int v : initialVertices){
            if(graph.getGraph().contains(v) || !triangleMap.containsKey(v))
                continue;

            // vertex v was deleted
            HashSet<Triangle> tri = triangleMap.remove(v);
            for(Triangle t : tri){
                for(int w : t.getVertices()){
                    if(triangleMap.containsKey(w))
                        triangleMap.get(w).remove(t);
                }
            }
        }
    }

    /**
     * Precomputed the sizes of the two-neighborhood of each vertex. Required for the order of start points.
     */
    private ArrayList<Integer> prepareStartPoints(){
        long time = System.currentTimeMillis();
        boolean aborted = false;
        startPointSizes = new int[this.highestVertexID + 1];
        ArrayList<Integer> ids = new ArrayList<>();
        for(int id: graph.getGraph().getVertices()){
            if(timeLimitEnabled && System.currentTimeMillis() > abortTime && ids.size() > 0){
                print("skip analyzing start points, time limit reached!");
                aborted = true;
                break;
            }
            startPointSizes[id] = graph.getGraph().sizeOfTwoNeighborhood(id, true);
            ids.add(id);
        }
        if(!aborted){
            ids = startPointQuickSort(ids);
            time = System.currentTimeMillis() - time;
        }
        output.addInformation("StartPointPreparationAborted", String.valueOf(aborted));
        output.setFormat("StartPointPreparationTime[s]", Output.DataFormat.MILLISECONDS);
        output.log("StartPointPreparation", time);
        return ids;
    }

    private void computeTriangleMap(){
        print("computing triangles ...");
        triangleMap = Triangle.getTriangleHashMap(graph.getGraph(), abortTime - System.currentTimeMillis());
        if(triangleMap == null)
            output.addInformation("TriangleComputationAborted", String.valueOf(true));
    }

    /**
     * Sorts the start points based on their two-neighborhood sizes in ascending order.
     * The sizes of their two neighborhoods are precomputed in the object array startPointSizes.
     * @param toSort the list of possible start points
     * @return the ordered list
     */
    private ArrayList<Integer> startPointQuickSort(ArrayList<Integer> toSort){
        if(toSort.isEmpty())
            return toSort;
        int pivotElem = toSort.get((int) (Math.random() * (toSort.size())));
        int pivotVal = startPointSizes[pivotElem];
        ArrayList<Integer> equalsList = new ArrayList<Integer>();
        ArrayList<Integer> smallerList = new ArrayList<Integer>();
        ArrayList<Integer> greaterList = new ArrayList<Integer>();
        for (int elem : toSort) {
            if (startPointSizes[elem] < pivotVal) smallerList.add(elem);
            else if (startPointSizes[elem] > pivotVal) greaterList.add(elem);
            else equalsList.add(elem);
        }
        ArrayList<Integer> solution = new ArrayList<Integer>();
        if (!smallerList.isEmpty()) solution.addAll(startPointQuickSort(smallerList));
        solution.addAll(equalsList);
        if (!greaterList.isEmpty()) solution.addAll(startPointQuickSort(greaterList));
        return solution;
    }

    /**
     * Sets formats and some graph information in the output.
     * @param graphName the name of the graph
     */
    private void prepareOutput(String graphName){
        output.addInformation("GraphName", graphName);
        output.addInformation("VertexCount", String.valueOf(inputGraph.size()));
        output.addInformation("EdgeCount", String.valueOf(inputGraph.getEdgeCount()));
        output.addInformation("TimeLimitEnabled", String.valueOf(timeLimitEnabled));
        output.addInformation("EdgeVariant", String.valueOf(useEdgeCondition));
        output.addInformation("L", String.valueOf(ell));
        output.addInformation("BranchingAborted", String.valueOf(false));
        output.addInformation("TriangleComputationAborted", String.valueOf(false));
        output.setFormat("Time[s]", Output.DataFormat.MILLISECONDS);
        output.setFormat("PreProcTime[s]", Output.DataFormat.MILLISECONDS);
        output.setFormat("StartPointCount", Output.DataFormat.INTEGER);
        output.setFormat("CurrentStartPoint", Output.DataFormat.INTEGER);
        output.setFormat("StartPointSkips", Output.DataFormat.INTEGER);
        output.setFormat("UpperBoundSkips", Output.DataFormat.INTEGER);
        output.setFormat("BranchCount", Output.DataFormat.INTEGER);
        output.setFormat("MaxBranchDepth", Output.DataFormat.INTEGER);
        output.setFormat("TimeLimit[s]", Output.DataFormat.MILLISECONDS);
        output.setFormat("SolutionSize", Output.DataFormat.INTEGER);

        output.log("TimeLimit[s]", algorithmTimeLimit);
        output.log("MaxBranchDepth", 0);
        output.log("StartPointCount", 0);
        output.log("StartPointSkips", 0);
        output.log("UpperBoundSkips", 0);
        output.log("BranchCount", 0);
    }

    /**
     * Prints something to the console if the print parameter in the run method was set
     * @param line a line to print
     */
    private void print(String line){
        if(printInformationToConsole)
            System.out.println(line);
    }

    private void resetAlgorithm(){
        output = new Output();
        conditionChecker = new ConditionChecker(this);
        graph = null;
        inputGraph = null;
        abortCurrentBranch = false;
        abortCurrentRun = false;
        timeLimitReached = false;
        triangleMap = null;
        startPointSizes = null;
        upperBound = null;
    }

    /**
     * Used to apply all data reduction rules before returning solution.
     */
    private static void applyDataReductionToGraph(RtcGraph g, boolean edgeCondition, int l){
        if(g == null)
            return;
        HashMap<Integer, HashSet<Triangle>> t = Triangle.getTriangleHashMap(g);
        Output o = new Output();
        ITriangleTwoCliqueAlgorithm a = new ITriangleTwoCliqueAlgorithm() {
            @Override
            public Map<Integer, Set<Integer>> identifyVerticesWithDistanceMaxTwo(T2CGraph graph) {
                return null;
            }

            @Override
            public Set<Set<Integer>> getSetsOfTriangle() {
                return null;
            }

            @Override
            public Map<Integer, Set<Set<Integer>>> mapRemainingTrianglesInGraph(T2CGraph globalGraph) {
                return null;
            }


            @Override
            public boolean containsOnlyTriangle(RtcGraph graph) {
                return false;
            }

            @Override
            public boolean isClique(RtcGraph graph) {
                return false;
            }

            @Override
            public Set getElementsOfSetOneNotInSetTwo(Set setOne, Set setTwo) {
                return null;
            }

            @Override
            public Set overlapOfTwoSets(Set setA, Set setB) {
                return null;
            }

            @Override
            public T2CGraph getOriginalGraph() {
                return null;
            }

            @Override
            public RtcGraph getGraph() {
                return g;
            }

            @Override
            public RtcGraph getConflictGraph() {
                return null;
            }

            @Override
            public HashSet<Integer> getMarked() {
                return null;
            }

            @Override
            public HashMap<Integer, HashSet<Triangle>> getTriangles() {
                return t;
            }

            @Override
            public ConditionChecker getConditionChecker() {
                return new ConditionChecker(this);
            }

            @Override
            public Output getOutput() {
                return o;
            }

            @Override
            public int getL() {
                return l;
            }

            @Override
            public int getBestSolutionSize() {
                return g.size();
            }

            @Override
            public boolean isEdgeConditionEnabled() {
                return edgeCondition;
            }

            @Override
            public void abortBranch() {}

            @Override
            public void deleteVertex(int v) {
                g.deleteVertex(v);
            }

            @Override
            public void deleteEdge(int v, int w) {
                g.getNeighbors(v).remove(w);
                g.getNeighbors(w).remove(v);
            }

            @Override
            public void markVertex(int v) {}
        };
        ArrayList<DeleteRule> rules = new ArrayList<>();
        rules.add(new VertexDegreeRule(a));
        if(edgeCondition)
            rules.add(new EdgeTriangleRule(a));
        rules.add(new VertexTriangleRule(a));

        boolean verticesDeleted;
        do{
            verticesDeleted = false;
            for(DataReductionRule d : rules){
                // skip other rules
                if(!(d instanceof DeleteRule)) continue;

                Double oldCount = o.getLoggedValue(d.getRuleName());
                d.apply();

                // check if vertices were deleted
                if(!o.getLoggedValue(d.getRuleName()).equals(oldCount))
                    verticesDeleted = true;

            }
        }while(verticesDeleted);
    }

    /*
     *  Implementation of interface methods
     */

    @Override
    public RtcGraph getGraph() {
        return graph == null ? null : graph.getGraph();
    }

    @Override
    public RtcGraph getConflictGraph() {
        return graph == null ? null : graph.getConflictGraph();
    }

    @Override
    public HashSet<Integer> getMarked() {
        return graph.getMarkedVertices();
    }

    @Override
    public HashMap<Integer, HashSet<Triangle>> getTriangles() {
        return triangleMap;
    }

    @Override
    public ConditionChecker getConditionChecker() {
        return this.conditionChecker;
    }

    @Override
    public Output getOutput() {
        return this.output;
    }

    @Override
    public int getL() {
        return ell;
    }

    @Override
    public int getBestSolutionSize() {
        return currentBestSolution == null ? 0 : currentBestSolution.size();
    }

    @Override
    public boolean isEdgeConditionEnabled(){
        return this.useEdgeCondition;
    }

    @Override
    public void abortBranch() {
        abortCurrentBranch = true;
    }

    @Override
    public void deleteVertex(int v) {
        if(graph.getGraph().contains(v))
            graph.deleteVertex(v);
    }

    @Override
    public void deleteEdge(int v, int w) {
        graph.deleteEdge(v, w);
    }

    @Override
    public void markVertex(int v) {
        graph.markVertex(v);
    }
}
