/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

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
     * Log levels
     */
    public static enum Level {
        /**
         * No logging
         */
        NONE(0),
        /**
         * Log only error messages
         */
        ERROR(100),
        /**
         * Log error and info messages
         */
        INFO(200),
        /**
         * Log debug messages
         */
        DEBUG(300);
        /**
         * Integer representation
         */
        final int level;

        Level(int level) {
            this.level = level;
        }
    }

    /**
     * the log level
     */
    private static Level logLevel = Level.INFO;

    /**
     * The console progress
     */
    private static final Progressable progress = new PrintstreamProgressable(System.err);

    /**
     * Is the logger interactive
     */
    private static boolean interactive = true;

    /**
     * Translates a string to the proper log level or throws an {@link IllegalArgumentException}
     *
     * @param level the current log level
     */
    public static void setLogLevel(String level) {
        try {
            logLevel = Level.valueOf(level.toUpperCase());
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Unknown log level " + level + ". Must be one of NONE|INFO|ERROR|DEBUG");
        }
    }

    /**
     * Set the log level
     *
     * @param level the level
     */
    public static void setLogLevel(Level level) {
        if (level == null) {
            throw new NullPointerException("Null LogLevel not permitted!");
        }
        logLevel = level;
    }

    /**
     * Get the current log level
     *
     * @return level the log level
     */
    public static Level getLogLevel() {
        return logLevel;
    }


    /**
     * Log a debug message
     *
     * @param message the message
     */
    public static void debug(String message) {
        debug(message, null);
    }

    /**
     * Log a debug message and print the stacktrace of the error
     *
     * @param message the message
     * @param error   the error (null permitted)
     */
    public static void debug(String message, Throwable error) {
        if (logLevel.level >= Log.Level.DEBUG.level) {
            System.err.println("[DEBUG] " + message);
            if (error != null) {
                error.printStackTrace(System.err);
            }
        }
    }

    /**
     * Log an waring message
     *
     * @param message the message
     */
    public static void warn(String message) {
        warn("WARN", message);
    }

    /**
     * Print a waring with a custom prefix
     *
     * @param prefix  the prefix
     * @param message the message
     */
    public static void warn(String prefix, String message) {
        info(prefix, message);
    }


    /**
     * Log an info message
     *
     * @param message the message
     */
    public static void info(String message) {
        info("INFO", message);
    }

    /**
     * Log an info message using a given prefix
     *
     * @param prefix prefix
     * @param message the message
     */
    public static void info(String prefix, String message) {
        if (logLevel.level >= Log.Level.INFO.level) {
            if (prefix != null && !prefix.isEmpty()) {
                System.err.println("[" + prefix + "] " + message);
            } else {
                System.err.println(message);
            }
        }
    }


    /**
     * Log an error message
     *
     * @param message the message
     */
    public static void error(String message) {
        error("ERROR", message);
    }

    /**
     * Log an error message with a custom error prefix
     *
     * @param prefix  the error prefix (null permitted)
     * @param message the message
     */
    public static void error(String prefix, String message) {
        if (logLevel.level >= Log.Level.ERROR.level) {
            if (prefix != null && prefix.length() > 0) {
                System.err.println("[" + prefix + "]" + " " + message);
            } else {
                System.err.println(message);
            }
        }
    }

    /**
     * Log an error message and print the stacktrace of the error
     *
     * @param message the message
     * @param error   the error
     */
    public static void error(String message, Throwable error) {
        if (logLevel.level >= Log.Level.ERROR.level) {
            error(message);
            if (error != null) {
                error.printStackTrace(System.err);
            }
        }
    }


    /**
     * Print the message as long as we are not in silent mode. The message is printed 'as is'
     * without any prefix.
     *
     * @param message the message
     */
    public static void message(String message) {
        if (logLevel.level >= Log.Level.INFO.level) {
            System.err.println(message);
        }
    }

    /**
     * Print to the default output directly. This does not check any log levels and always prints
     *
     * @param message the message
     */
    public static void println(String message) {
        System.out.println(message);
    }

    /**
     * Print to the default output directly. This does not check any log levels and always prints
     *
     * @param message the message
     */
    public static void print(String message) {
        System.out.print(message);
    }


    /**
     * Start a new progress
     *
     * @param message the message (null permitted)
     */
    public static void progressStart(String message) {
        if (logLevel.level >= Log.Level.INFO.level) {
            progress.start(message);
        }
    }

    /**
     * Do a progress step
     */
    private static void progress() {
        if (logLevel.level >= Log.Level.INFO.level) {
            progress.progress();
        }
    }

    /**
     * Checks if the progress should be printed
     *
     * @param currentValue the current value
     * @param maxValue     the maximum value
     */
    public static void progress(long currentValue, long maxValue) {
        if (logLevel.level >= Log.Level.INFO.level) {
            int preStep = (int) Math.ceil((double) maxValue / progress.steps());
            while (currentValue < maxValue && progress.currentStep() < progress.steps() && progress.currentStep() * preStep <= currentValue) {
                progress();
            }
        }
    }


    /**
     * Finish the current progress
     */
    public static void progressFinish() {
        if (logLevel.level >= Log.Level.INFO.level) {
            progress.finish();
        }
    }

    /**
     * Finish the current progress with an optional message and optionally print the time
     *
     * @param msg  the message (null permitted)
     * @param time print the time
     */
    public static void progressFinish(String msg, boolean time) {
        if (logLevel.level >= Log.Level.INFO.level) {
            progress.finish(msg, time);
        }
    }

    /**
     * Progress failed
     *
     * @param msg the message (null permitted)
     */
    public static void progressFailed(String msg) {
        if (logLevel.level >= Log.Level.INFO.level) {
            progress.failed(msg);
        }
    }

    /**
     * Returns the number of steps provided by the progress
     *
     * @return steps the number of steps provided by the progress
     */
    public static int progressSteps() {
        return progress.steps();
    }

    /**
     * Returns true if this logger is interactive and user can input things
     *
     * @return interactive interactive
     */
    public static boolean isInteractive() {
        return interactive;
    }

    /**
     * Set this logger interactive. If true, the user can input things
     *
     * @param interactive interactive
     */
    public static void setInteractive(boolean interactive) {
        Log.interactive = interactive;
    }
}
