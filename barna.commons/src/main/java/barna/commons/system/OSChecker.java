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

package barna.commons.system;

/**
 * Check the current operating system
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class OSChecker {

    private static boolean win7= false;
    private static boolean winVista= false;
    private static boolean winXP= false;
    private static boolean win2K= false;
    private static boolean winME= false;
    private static boolean winNT= false;
    private static boolean win9x= false;
    private static boolean linux= false;
    private static boolean sunOS= false;
    private static boolean hpOS= false;
    private static boolean macOSX= false;
    private static boolean is32bit= false;
    private static boolean is64bit= false;
    /**
     * The new line separator
     */
    public static String NEW_LINE = System.getProperty("line.separator");

    public static enum OS {
        win("Windows"), linux("Linux"), sun("Solaris (x86)"), mac("Mac OS"), unknown("Unknown");

        private String name;

        OS(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }
    }

    static {
        init();
    }


    private static void init() {

        String osName= System.getProperty("os.name").toLowerCase();
        if (osName.indexOf("windows 7") > -1)
            win7= true;
        if (osName.indexOf("windows vista") > -1)
            winVista = true;
        if (osName.indexOf("windows xp") > -1)
            winXP= true;
        if (osName.indexOf("windows 2000") > -1)
            win2K= true;
        if (osName.indexOf("windows me") > -1)
            winME= true;
        if (osName.indexOf("windows nt") > -1)
            winNT= true;
        if (osName.indexOf("windows 9") > -1)
            win9x= true;

        if (osName.indexOf("linux") > -1) 		// startsWith()
            linux= true;
        if ((osName.indexOf("sunos") > -1)||
            (osName.indexOf("solaris") > -1)) 	// startsWith()
            sunOS= true;
        if (osName.indexOf("hp") > -1) 		// startsWith()
            hpOS= true;

        if (osName.indexOf("mac os x") > -1) 		// startsWith()
            macOSX= true;

        String osArch = System.getProperty("os.arch");
        if (osArch.indexOf("64") > -1) {
            is64bit = true;
        } else {
            is32bit = true;
        }
    }


    /**
     * Returns the linux.
     * @return boolean
     */
    public static boolean isLinux() {
        return linux;
    }

    /**
     * Returns the sunOS.
     * @return boolean
     */
    public static boolean isSunOS() {
        return sunOS;
    }

    /**
     * Returns the win2K.
     * @return boolean
     */
    public static boolean isWin2K() {
        return win2K;
    }

    /**
     * Returns the win9x.
     * @return boolean
     */
    public static boolean isWin9x() {
        return win9x;
    }

    /**
     * Returns the winME.
     * @return boolean
     */
    public static boolean isWinME() {
        return winME;
    }

    /**
     * Returns the winNT.
     * @return boolean
     */
    public static boolean isWinNT() {
        return winNT;
    }

    /**
     * Returns the winXP.
     * @return boolean
     */
    public static boolean isWinXP() {
        return winXP;
    }

    /**
     * Returns the winVista
     *
     * @return boolean
     */
    public static boolean isWinVista() {
        return winVista;
    }

    /**
     * Returns the win7
     * @return boolean
     */
    public static boolean isWin7() {
        return win7;
    }

    /**
     * Returns the windows.
     * @return boolean
     */
    public static boolean isWindows() {
        return (isWinXP()|| isWin2K()|| isWin9x()|| isWinME()|| isWinNT() || isWinVista() || isWin7());
    }

    /**
     * Returns the macOSX.
     * @return boolean
     */
    public static boolean isMacOSX() {
        return macOSX;
    }

    /**
     * Returns the hpOS.
     * @return boolean
     */
    public static boolean isHpOS() {
        return hpOS;
    }
    /**
     * Get's the version of Java currently running.
     *
     * @return the version of Java that is running.
     */
    public static String getJavaVersion() {
        return System.getProperty("java.version");
    }

    /**
     * Gets the operating system version that the JVM is running on.
     *
     * @return the operating system version that the JVM is running on.
     */
    public static String getOsVersion() {
        return System.getProperty("os.version");
    }

    /**
     * True if this JVM is running on a Mac.
     *
     * @return true if this JVM is running on a Mac.
     */
    public static boolean isMac() {
        return System.getProperty("os.name").startsWith("Mac OS");
    }

    /**
     * True if this JVM is running Java 6 on a Mac.
     *
     * @return true if this JVM is running Java 6 on a Mac.
     */
    public static boolean isJava6OnMac() {
        return isMac() && getJavaVersion().startsWith("1.6");
    }

    /**
     * True if this JVM is running 64 bit Java on a Mac.
     *
     * @return true if this JVM is running 64 bit Java on a Mac.
     */
    public static boolean is64BitJavaOnMac() {
        return isMac() && System.getProperty("os.arch").equals("x86_64");
    }

    /**
     * True if this JVM is running on Mac OS X 10.5, Leopard.
     *
     * @return true if this JVM is running on Mac OS X 10.5, Leopard.
     */
    public static boolean isLeopard() {
        return isMacOSX() && getOsVersion().startsWith("10.5");
    }

    /**
     * True if this JVM is running 32 bit Java.
     *
     * @return true if this JVM is running 32 bit Java.
     */
    public static boolean is32bit() {
        return is32bit;
    }

    /**
     * True if this JVM is running 64 bit Java.
     *
     * @return true if this JVM is running 64 bit Java.
     */
    public static boolean is64bit() {
        return is64bit;
    }

    /**
     * Get the current OS
     *
     * @return os the current os
     */
    public static OS getOS(){
        if(isWindows()) return OS.win;
        if(isLinux()) return  OS.linux;
        if(isMac() || isMacOSX()) return OS.mac;
        if(isSunOS()) return OS.sun;
        return OS.unknown;
    }

}
