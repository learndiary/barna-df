package fbi.commons;

/**
 * Step based progress. Use {@link #steps()}  to find the number of steps provided
 * by the implementation. Call {@link #start(String)} to reset the current steps and start
 * a new progress. If a current progress is running, it is canceled without a message.
 */
public interface Progressable {
    /**
     * Start a new progress run. This ends running progresses and resets the step counter.
     * Implementation must also store the start time to be able to print the start time
     * when the progress is finished
     *
     * @param message the message (null permitted)
     */
    public void start(String message);

    /**
     * Do one progress step if more steps are available
     */
	public void progress();

    /**
     * Finish the progress. This print the rest of the available steps and then ends the progress
     */
	public void finish();

    /**
     * Finish the progress and print the optional message
     *
     * @param msg the message (null permitted)
     * @param time print the running time
     */
	public void finish(String msg, boolean time);

    /**
     * Progress process failed. This ends the progress and prints the (optional) message
     *
     * @param msg the message (null permitted)
     */
    public void failed(String msg);

    /**
     * Returns the number of steps provided by this progressable.
     *
     * @return steps the number of steps
     */
    public int steps();

    /**
     * Returns the current step
     *
     * @return step the current step
     */
    public int currentStep();
}
