package fbi.commons.tools;

import fbi.commons.Log;

/**
 * Static helper methods for command line tools
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class CommandLine {
    /**
     * Yes options for user confirmation
     */
    public static final String[] user_yes = new String[]{"yes", "y", "si", "yo", "ja", "ok"};
    /**
     * No options for user confirmation
     */
    public static final String[] user_no = new String[]{"no", "n", "nein", "nope", "nix", "noe"};


    /**
     * Ask the user to confirm the given question. If logging is disabled (shutup mode), this will always return true
     *
     * @param message the question (null permitted)
     * @return confirmed true if the user confirmed with yes
     */
    public static boolean confirm(String message) {
        if (Log.getLogLevel() == Log.Level.NONE || !Log.isInteractive())
            return true;

        while (true) {
            StringBuffer sb = new StringBuffer();
            int in;

            if (message != null) {
                System.err.print(message + " ");
                System.err.flush();
            }

            try {
                while ((in = System.in.read()) != '\n')
                    sb.append((char) in);
            } catch (Exception e) {
                ; // :)
            }

            String s = sb.toString().toLowerCase().trim();
            for (int i = 0; i < user_yes.length; i++)
                if (s.equals(user_yes[i]))
                    return true;
            for (int i = 0; i < user_no.length; i++)
                if (s.equals(user_no[i]))
                    return false;
        }
    }

}
