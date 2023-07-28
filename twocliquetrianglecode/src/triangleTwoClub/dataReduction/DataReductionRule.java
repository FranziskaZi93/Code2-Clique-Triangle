package triangleTwoClub.dataReduction;

import triangleTwoClub.Output;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

public abstract class DataReductionRule {

    ITriangleTwoClubAlgorithm algorithm;
    private boolean isRegistered = false;
    private boolean isActivated = true;
    private String timeLogName;
    private String callCountName;

    public DataReductionRule(ITriangleTwoClubAlgorithm algorithm){
        this.algorithm = algorithm;
    }

    public void apply(){
        if(!isRegistered)
            registerOutput();
        if(!isActivated)
            return;

        algorithm.getOutput().log(callCountName, 1L);
        long time = System.currentTimeMillis();
        applyRule();
        algorithm.getOutput().log(timeLogName, System.currentTimeMillis() - time);
    }

    private void registerOutput(){
        // statistic
        algorithm.getOutput().setFormat(getRuleName(), Output.DataFormat.INTEGER);
        algorithm.getOutput().log(getRuleName(), 0);

        // time measurement
        timeLogName = getRuleName() + "Time[s]";
        algorithm.getOutput().setFormat(timeLogName, Output.DataFormat.MILLISECONDS);

        // call count
        callCountName = getRuleName() + "Count";
        algorithm.getOutput().setFormat(callCountName, Output.DataFormat.INTEGER);
        algorithm.getOutput().log(callCountName, 0);

        isRegistered = true;
    }

    public void setActivated(boolean useRule){
        this.isActivated = useRule;
    }

    public abstract String getRuleName();
    protected abstract void applyRule();
    public abstract boolean applyInitially();
}
