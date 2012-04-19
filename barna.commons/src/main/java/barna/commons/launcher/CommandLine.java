/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.commons.launcher;

import barna.commons.log.Log;

/**
 * Static helper methods for command line tools. This class is bound to the {@link Log} implementation
 * to check for quiet state and log levels.
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
        if (Log.getLogLevel() == Log.Level.NONE || !Log.isInteractive()) {
            return true;
        }

        while (true) {
            StringBuilder sb = new StringBuilder();
            int in;

            if (message != null) {
                System.err.print(message + " ");
                System.err.flush();
            }

            try {
                while ((in = System.in.read()) != '\n') {
                    sb.append((char) in);
                }
            } catch (Exception ignore) {
                // ignore exceptions
            }

            String s = sb.toString().toLowerCase().trim();
            for (String user_ye : user_yes) {
                if (s.equals(user_ye)) {
                    return true;
                }
            }
            for (String anUser_no : user_no) {
                if (s.equals(anUser_no)) {
                    return false;
                }
            }
        }
    }

}
