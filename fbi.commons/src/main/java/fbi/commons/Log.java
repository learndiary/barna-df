package fbi.commons;

/**
 * Defines constants and current logging level.
 * This is a simple abstraction for the system out stream at the moment, but we
 * might incorporate a more powerful logging mechanism later.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Log {
    /**
     * Debugging mode
     */
    public static final byte VERBOSE_DEBUG = 3;
    /**
     * Default log level
     */
    public static final byte VERBOSE_NORMAL = 2;
    /**
     * Show errors only
     */
    public static final byte VERBOSE_ERRORS = 1;
    /**
     * Quite log level
     */
    public static final byte VERBOSE_SHUTUP = 0; // todo: rename to "quiet"
    /**
     * String representations for logging levels
     */
    public static final String[] VERBOSE_KEYWORDS= new String[] {"SILENT", "VERBOSE", "ERRORS", "DEBUG"};
    /**
     * Current logging level
     */
    public static byte verboseLevel= VERBOSE_NORMAL; // todo: does this need to be public ?

    /**
     * The console progress
     */
    private static Progressable progress = new PrintstreamProgressable(System.err);

    /**
     * Set logging level and returns true if the given string was a proper logging level, otherwise this returns false
     * and the logging level is not changed
     *
     * @param logLevel the current log level
     * @return changed returns true if log level was applied successfully
     */
    public static final boolean setVerbose(String logLevel) {
        //SILENT, VERBOSE, ERRORS, DEBUG
        logLevel= logLevel.toUpperCase();
        for (int i = 0; i < VERBOSE_KEYWORDS.length; i++){
            if (logLevel.equals(VERBOSE_KEYWORDS[i])) {
                verboseLevel= (byte) i;
                return true;
            }
        }
        return false;
    }


    /**
     * Log a debug message
     *
     * @param message the message
     */
    public static void debug(String message){
        if(verboseLevel >= VERBOSE_DEBUG)
            System.err.println("[DEBUG] " + message);
    }

    /**
     * Log a debug message and print the stacktrace of the error
     *
     * @param message the message
     * @param error the error
     */
    public static void debug(String message, Throwable error){
        if(verboseLevel >= VERBOSE_DEBUG){
            debug(message);
            if(error != null)error.printStackTrace(System.err);
        }
    }


    /**
     * Log an info message
     *
     * @param message the message
     */
    public static void info(String message){
        info("INFO", message);
    }

    /**
     * Log an waring message
     *
     * @param message the message
     */
    public static void warn(String message){
        info("WARN", message);
    }

    public static void warn(String prefix, String message){
        info("WARN", message);
    }



    /**
     * Log an info message using a given prefix
     *
     * @param message the message
     */
    public static void info(String prefix, String message){
        if(verboseLevel >= VERBOSE_NORMAL){
            if(prefix != null && !prefix.isEmpty()){
                System.err.println("["+prefix+"] "+message);
            }else{
                System.err.println(message);
            }
        }
    }


    /**
     * Log an error message
     *
     * @param message the message
     */
    public static void error(String message){
        if(verboseLevel >= VERBOSE_ERRORS){
            System.err.println("[ERROR] "+message);
        }
    }
    /**
     * Log an error message and print the stacktrace of the error
     *
     * @param message the message
     * @param error the error
     */
    public static void error(String message, Throwable error){
        if(verboseLevel >= VERBOSE_ERRORS){
            error(message);
            if(error != null)error.printStackTrace(System.err);
        }
    }


    /**
     * Print the message as long as we are not in silent mode. The message is printed 'as is'
     * without any prefix.
     *
     * @param message the message
     */
    public static void message(String message){
        if(verboseLevel > VERBOSE_SHUTUP){
            System.err.println(message);
        }
    }

    /**
     * Start a new progress
     *
     * @param message the message (null permitted)
     */
    public static void progressStart(String message) {
        progress.start(message);
    }

    /**
     * Do a progress step
     */
    public static void progress(){
        progress.progress();
    }
    /**
     * Checks if the progress should be printed
     *
     * @param currentValue the current value
     * @param maxValue the maximum value
     */
    public static void progress(long currentValue, long maxValue) {
        int preStep = (int) Math.ceil((double) maxValue / progress.steps());
        while(currentValue < maxValue && progress.currentStep() < progress.steps() && progress.currentStep() * preStep <= currentValue){
            progress();
        }
    }


    /**
     * Finish the current progress
     */
    public static void progressFinish() {
        progress.finish();
    }

    /**
     * Finish the current progress with an optional message and optionally print the time
     *
     * @param msg the message (null permitted)
     * @param time print the time
     */
    public static void progressFinish(String msg, boolean time) {
        progress.finish(msg, time);
    }

    /**
     * Progress failed
     *
     * @param msg the message (null permitted)
     */
    public static void progressFailed(String msg) {
        progress.failed(msg);
    }

    /**
     * Returns the number of steps provided by the progress
     *
     * @return steps the number of steps provided by the progress
     */
    public static int progressSteps() {
        return progress.steps();
    }
}
