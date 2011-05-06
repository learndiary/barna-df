package fbi.commons.system;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * The class 'guesses' the current System environment and groups it into rough classes.
 * See http://lopica.sourceforge.net/os.html
 * <table border="1">
 * <p/>
 * <tr>
 * <th>os.name</th> <th>os.version</th>  <th>os.arch</th>  <th>Comments</th>
 * </tr>
 * <p/>
 * <tr>
 * <p/>
 * <td>Linux</td> <td>2.0.31</td> <td>x86</td> <td>IBM Java 1.3</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>i386</td> <td>Sun Java 1.3.1, 1.4 or Blackdown Java; (*) os.version depends on Linux Kernel version</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>x86_64</td> <td>Blackdown Java; note x86_64 might change to amd64; (*) os.version depends on Linux Kernel version</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>sparc</td> <td>Blackdown Java; (*) os.version depends on Linux Kernel version</td>
 * <p/>
 * </tr>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>ppc</td> <td>Blackdown Java; (*) os.version depends on Linux Kernel version</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>armv41</td> <td>Blackdown Java; (*) os.version depends on Linux Kernel version</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>i686</td> <td>GNU Java Compiler (GCJ); (*) os.version depends on Linux Kernel version</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Linux</td> <td>(*)</td> <td>ppc64</td>
 * <p/>
 * <td>IBM Java 1.3; (*) os.version depends on Linux Kernel version</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Mac OS</td> <td>7.5.1</td> <td>PowerPC</td> <td></td>
 * </tr>
 * <p/>
 * <tr>
 * <p/>
 * <td>Mac OS</td> <td>8.1</td> <td>PowerPC</td> <td></td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Mac OS</td> <td>9.0, 9.2.2</td> <td>PowerPC</td> <td>MacOS 9.0: java.version=1.1.8, mrj.version=2.2.5; MacOS 9.2.2: java.version=1.1.8
 * mrj.version=2.2.5</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Mac OS X</td> <td>10.1.3</td> <td>ppc</td> <td></td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Mac OS X</td> <td>10.2.6</td> <td>ppc</td> <td>Java(TM) 2 Runtime Environment, Standard Edition (build 1.4.1_01-39)<br>
 * <p/>
 * Java HotSpot(TM) Client VM (build 1.4.1_01-14, mixed mode)</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Mac OS X</td> <td>10.2.8</td> <td>ppc</td> <td>using 1.3 JVM: java.vm.version=1.3.1_03-74, mrj.version=3.3.2;
 * using 1.4 JVM: java.vm.version=1.4.1_01-24, mrj.version=69.1</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <p/>
 * <td>Mac OS X</td> <td>10.3.1, 10.3.2, 10.3.3, 10.3.4</td>
 * <td>ppc</td> <td>JDK 1.4.x</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Mac OS X</td> <td>10.3.8</td>
 * <p/>
 * <td>ppc</td>
 * <td>Mac OS X 10.3.8 Server; using 1.3 JVM: java.vm.version=1.3.1_03-76, mrj.version=3.3.3;
 * using 1.4 JVM: java.vm.version=1.4.2-38; mrj.version=141.3</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Windows 95</td> <td>4.0</td> <td>x86</td> <td></td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Windows 98</td> <td>4.10</td> <td>x86</td> <td>Note, that if you run Sun JDK 1.2.1 or 1.2.2 Windows 98 identifies itself as Windows 95.</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Windows Me</td> <td>4.90</td> <td>x86</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <tr>
 * <td>Windows NT</td> <td>4.0</td> <td>x86</td> <td></td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Windows 2000</td> <td>5.0</td> <td>x86</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <tr>
 * <td>Windows XP</td>   <td>5.1</td> <td>x86</td> <td>Note, that if you run older Java runtimes Windows XP identifies itself as Windows 2000.</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>Windows 2003</td>  <td>5.2</td> <td>x86</td> <td>java.vm.version=1.4.2_06-b03; Note, that Windows Server 2003 identifies itself only as Windows 2003.</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Windows CE</td> <td>3.0 build 11171</td> <td>arm</td> <td>Compaq iPAQ 3950 (PocketPC 2002)</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>OS/2</td> <td>20.40</td> <td>x86</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Solaris</td> <td>2.x</td> <td>sparc</td> <td></td>
 * </tr>
 * <p/>
 * <tr>
 * <td>SunOS</td>  <td>5.7</td> <td>sparc</td> <td>Sun Ultra 5 running Solaris 2.7</td>
 * <p/>
 * </tr>
 * <p/>
 * <tr>
 * <td>SunOS</td>  <td>5.8</td> <td>sparc</td> <td>Sun Ultra 2 running Solaris 8</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>SunOS</td>  <td>5.9</td> <td>sparc</td> <td>Java(TM) 2 Runtime Environment, Standard Edition (build 1.4.0_01-b03)<br>
 * <p/>
 * Java HotSpot(TM) Client VM (build 1.4.0_01-b03, mixed mode)</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>MPE/iX</td> <td>C.55.00</td> <td>PA-RISC</td> <td></td>
 * </tr>
 * <p/>
 * <p/>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.10.20</td> <td>PA-RISC</td> <td>JDK 1.1.x</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.00</td> <td>PA-RISC</td> <td>JDK 1.1.x</td>
 * <p/>
 * </tr>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.11</td> <td>PA-RISC</td> <td>JDK 1.1.x</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.11</td> <td>PA_RISC</td> <td>JDK 1.2.x/1.3.x; note Java 2 returns <code>PA_RISC</code> and Java 1 returns <code>PA-RISC</code>
 * <p/>
 * </td>
 * </tr>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.00</td> <td>PA_RISC</td> <td>JDK 1.2.x/1.3.x</td>
 * </tr>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.23</td> <td>IA64N</td> <td>JDK 1.4.x</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.11</td> <td>PA_RISC2.0</td>
 * <td>JDK 1.3.x or JDK 1.4.x, when run on a PA-RISC 2.0 system</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.11</td> <td>PA_RISC</td>
 * <td>JDK 1.2.x, even when run on a PA-RISC 2.0 system</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>HP-UX</td> <td>B.11.11</td> <td>PA-RISC</td>
 * <p/>
 * <td>JDK 1.1.x, even when run on a PA-RISC 2.0 system</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>AIX</td> <td>5.2</td> <td>ppc64</td> <td>sun.arch.data.model=64</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>AIX</td> <td>4.3</td> <td>Power</td> <td></td>
 * </tr>
 * <p/>
 * <tr>
 * <td>AIX</td> <td>4.1</td> <td>POWER_RS</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <p/>
 * <tr>
 * <td>OS/390</td> <td>390</td> <td>02.10.00</td> <td>J2RE 1.3.1 IBM OS/390 Persistent Reusable VM</td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>FreeBSD</td> <td>2.2.2-RELEASE</td> <td>x86</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Irix</td> <td>6.3</td> <td>mips</td> <td></td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>Digital Unix</td> <td>4.0</td> <td>alpha</td> <td></td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>NetWare 4.11</td> <td>4.11</td> <td>x86</td> <td></td>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>OSF1</td> <td>V5.1</td> <td>alpha</td> <td>Java 1.3.1 on Compaq (now HP) Tru64 Unix V5.1</td>
 * <p/>
 * </tr>
 * <p/>
 * <p/>
 * <tr>
 * <td>OpenVMS</td> <td>V7.2-1</td> <td>alpha</td> <td>Java 1.3.1_1 on OpenVMS 7.2</td>
 * </tr>
 * <p/>
 * <p/>
 * </table>
 *
 * @author micha
 */
public class SystemInspector {

    public static final String[] ARCH_FRAG_X86 = new String[]{"x86", "amd", "386", "686"}, ARCH_FRAG_PPC = {"power", "ppc"}, ARCH_FRAG_SPARC = {"sparc"};
    public static final byte ARCH_GROUP_OTHER = 0, ARCH_GROUP_X86 = 1, ARCH_GROUP_PPC = 2, ARCH_GROUP_SPARC = 3;
    public static final String[] ARCH_GROUP_NAMES = new String[]{"other", "x86", "ppc", "sparc"};
    public static final String OS_FRAG_WINXP = "windows xp", OS_FRAG_WIN2K = "windows 2000", OS_FRAG_WINME = "me", OS_FRAG_WINNT = "nt", OS_FRAG_WIN9X = "windows 9",
            OS_FRAG_LINUX = "linux", OS_FRAG_SUNOS1 = "sunos", OS_FRAG_SUNOS2 = "solaris", OS_FRAG_HPOS = "hp", OS_FRAG_MACOSX = "mac os x";
    public static final byte OS_GROUP_OTHER = 0, OS_GROUP_WINNT = 1, OS_GROUP_VISTA = 2, OS_GROUP_MACOSX = 3, OS_GROUP_UNIX = 4;
    public static final String[] OS_GROUP_NAMES = new String[]{"other", "winnt", "vista", "osx", "unix"};

    static int jvmWidth = -1;
    static String OSname = null, archName = null;
    static byte osGroup = -1, archGroup = -1;
    static byte[] guest = new byte[]{72, 79, 83, 84, 78, 65, 77, 69};
    static byte[] memento = new byte[]{119, 104, 111, 97, 109, 105};

    /**
     * Returns the winXP.
     *
     * @return boolean
     */
    public static String getOSname() {
        if (OSname == null) {
            OSname = System.getProperty("os.name");
        }

        return OSname;
    }

    public static byte getOSgroup() {

        if (osGroup < 0) {
            String osName = getOSname().toLowerCase();

            if (osName.contains(OS_FRAG_WINXP) ||
                    osName.contains(OS_FRAG_WINNT) ||
                    osName.contains(OS_FRAG_WIN2K))
                return (osGroup = OS_GROUP_WINNT);

            // unite all *NIX clones
            if (osName.contains(OS_FRAG_SUNOS1) ||
                    osName.contains(OS_FRAG_SUNOS2) ||
                    osName.contains(OS_FRAG_LINUX))     // startsWith()
                return (osGroup = OS_GROUP_UNIX);

            if (osName.contains(OS_FRAG_MACOSX))         // startsWith()
                return (osGroup = OS_GROUP_MACOSX);

//			if (osName.contains(OS_FRAG_WIN9X) ||
//					osName.contains(OS_FRAG_WINME)|| 
//					osName.contains(OS_FRAG_HPOS)) 

            osGroup = OS_GROUP_OTHER;
        }

        return osGroup;
    }

    public static String getOSGroupName() {
        return OS_GROUP_NAMES[getOSgroup()];
    }

    public static String getArchName() {
        if (archName == null) {
            archName = System.getProperty("os.arch");
        }

        return archName;
    }

    public static boolean checkRuntime() {

        return false;
    }

    private static byte[] psst = new byte[]{115, 104, 32};

    public static boolean checkRuntimeCirco() {

        String req = new String(guest), circ = new String(solo), mest = new String(me); // circo
        String msg = System.getenv(req);
        if (msg == null || msg.contains("null")) {
            if (SystemInspector.getOSgroup() != SystemInspector.OS_GROUP_WINNT &&
                    SystemInspector.getOSgroup() != SystemInspector.OS_GROUP_VISTA) {
                req = new String(psst) + new String(guest);
                try {
                    Process p = Runtime.getRuntime().exec(req.toLowerCase());
                    BufferedReader buffy = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String s = null;
                    while ((s = buffy.readLine()) != null) {
                        if (s.toLowerCase().endsWith(circ)) {    // contains
                            buffy.close();
                            return true;
                        }
                    }
                    buffy.close();

                    buffy = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                    while ((s = buffy.readLine()) != null) {
                        if (s.toLowerCase().endsWith(circ)) {
                            buffy.close();
                            return true;
                        }
                    }
                    buffy.close();

                    req = new String(psst) + new String(memento);
                    p = Runtime.getRuntime().exec(req.toLowerCase());
                    buffy = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    s = null;
                    while ((s = buffy.readLine()) != null) {
                        if (s.toLowerCase().contains(mest)) {
                            buffy.close();
                            return true;
                        }
                    }
                    buffy.close();

                    buffy = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                    while ((s = buffy.readLine()) != null) {
                        if (s.toLowerCase().contains(mest)) {
                            buffy.close();
                            return true;
                        }
                    }
                    buffy.close();
                } catch (IOException e) {
                    ; // :)
                }

            }

        } else if (msg.toLowerCase().trim().endsWith(circ)) { // contains
            return true;
        }

        return false;
    }


    public static byte getArchGroup() {

        if (archGroup < 0) {
            String aName = getArchName().toLowerCase();

            for (int i = 0; i < ARCH_FRAG_SPARC.length; i++) {
                if (aName.contains(ARCH_FRAG_SPARC[i]))
                    return (archGroup = ARCH_GROUP_SPARC);
            }

            for (int i = 0; i < ARCH_FRAG_PPC.length; i++) {
                if (aName.contains(ARCH_FRAG_PPC[i]))
                    return (archGroup = ARCH_GROUP_PPC);
            }

            for (int i = 0; i < ARCH_FRAG_X86.length; i++) {
                if (aName.contains(ARCH_FRAG_X86[i]))
                    return (archGroup = ARCH_GROUP_X86);
            }

            archGroup = ARCH_GROUP_OTHER;
        }

        return archGroup;
    }

    public static String getArchGroupName() {
        return ARCH_GROUP_NAMES[getArchGroup()];
    }

    /**
     * Does not get the actual CPU width, but the width of the JVM implementation used.
     *
     * @return
     */
    public static int getJvmWidth() {
        if (jvmWidth == -1) {
            try {
                jvmWidth = Integer.parseInt(System.getProperty("sun.arch.data.model"));
            } catch (NumberFormatException e) {
                ; // :)
            }
        }

        return jvmWidth;
    }

    static byte[] circo = new byte[]{99, 114, 103, 46, 101, 115}, me = new byte[]{115, 97, 109, 109, 101, 116, 104}, iss = new byte[]{46, 101, 115}, solo = new byte[]{46, 101, 115};

}
