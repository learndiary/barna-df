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

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.log.Log;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.io.File;
import java.io.FilePermission;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.security.AccessControlException;
import java.security.Permission;
import java.util.*;


/**
 * Flux Simulator starter class. Contains the main method and parses command line arguments. During startup,
 * this checks for any {@link Tool} implementations and adds them to the list of available utils.
 */
public class Flux {
    /**
     * Required java version
     */
    private static final float FLUX_JAVA_VERSION = 1.6f;
    /**
     * The tool to be executed
     */
    private String toolName;

    /**
     * Number of executor threads
     */
    private int threads = 2;



    /**
     * Start the Flux simulator
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
    	
    	// set a system manager
    	SecurityManager sm= getSecurityManager();
    	if (sm!= null)
    		System.setSecurityManager(sm);
    	
        // load java.util.logger configuration
        Log.initialize();

        // print tool and version info
        FluxVersionInfo versionInfo = null;
        if(System.getProperty("flux.app") != null){
            versionInfo = new FluxVersionInfo("/" + System.getProperty("flux.app")+"-build.properties");
            Log.println(versionInfo.toShortString());
        }


        /*
       Check java version
        */
        checkJavaVersion();
        JSAP jsap = createBaseOptions();


        // register tools
        List<Tool> tools = findTools();

        // prepare the fluxInstance
        Flux fluxInstance = new Flux();

        // parse the arguments (first round)
        JSAPResult initialFluxArguments = null;
        try{
            initialFluxArguments = jsap.parse(args);
        }catch(Exception e){
            Log.error("Error while parsing arguments : " + e.getMessage());
        }

        fluxInstance.setLogLevel(initialFluxArguments.getString("log"));
        fluxInstance.setThreads(initialFluxArguments.getInt("threads"));
        fluxInstance.setToolName(initialFluxArguments.getString("tool"));
        fluxInstance.setDetached(initialFluxArguments.userSpecified("force"));

        if(initialFluxArguments.userSpecified("version") && versionInfo != null){
            System.err.println(versionInfo.toString());
            System.exit(1);
        }

        // still no tool ? print usage and exit
        if(fluxInstance.getToolName() == null){
            printUsage(null, jsap, tools, null, false);
        }

        /**
         * Print tools
         */
        if(initialFluxArguments.userSpecified("list-tools")){
            printTools(tools);
            System.exit(-1);
        }


        // find the tool to start
        Tool tool = null;
        for (Tool fluxTool : tools) {
            if(fluxTool.getName().equals(fluxInstance.getToolName())){
                tool = fluxTool;
                break;
            }
        }

        /// register tool parameter
        if(tool != null){
            List<Parameter> parameter = tool.getParameter();
            if(parameter != null){
                try{
                    for (Parameter p : parameter) {
                        jsap.registerParameter(p);
                    }
                } catch (Exception e) {
                    Log.error("Parameter error : " + e.getMessage(), e);
                    System.exit(-1);
                }
            }
        }

        boolean userRequestsHelp = initialFluxArguments.userSpecified("help");
        if (userRequestsHelp || tool == null) {
            // todo: add error message "No tool sepcified"
            printUsage(tool, jsap, tools, !userRequestsHelp ?
                    (initialFluxArguments.userSpecified("tool") ?
                            initialFluxArguments.getString("tool") + " tool not found!" :
                            "No tool specified, use -t <tool> to specify a tool")
                    :null, userRequestsHelp);
        }

        // execute the tool
        // 1. get the parameters
        // 2. parse them
        // 3. let the tool validate the parameter
        List<Parameter> parameter = tool.getParameter();
        if(parameter != null){
            try{
                JSAPResult toolParameter = jsap.parse(args);
                if(!tool.validateParameter(toolParameter)){
                    printUsage(tool, jsap, tools, null, false);
                }
            } catch (Exception e) {
                Log.error("Parameter error : " + e.getMessage(), e);
                System.exit(-1);
            }

        }

        try {
            // configure the executor
            Execute.initialize(fluxInstance.getThreads());

            tool.call();
        }catch (OutOfMemoryError outOfMemoryError){
            int mb = 1024*1024;
            long maxMemoryBytes = Runtime.getRuntime().maxMemory();
            long maxMemoryMB = (maxMemoryBytes/mb);
            Log.error("Execution failed!\n"+
                    "The Flux " + fluxInstance.getToolName() + " tool run into memory problems !" +
                    "Please use the FLUX_MEM environment variable to increase the memory. For example: export FLUX_MEM=\"6G\"; " +
                    "flux-"+fluxInstance.getToolName()+" ... to use " +
                    "6 GB of memory.");
            Log.error("Current memory setting : " + maxMemoryMB + " MB");
            Log.error(("Tool that run into memory issues: " + tool.getName()));
            if(outOfMemoryError.getStackTrace() != null && outOfMemoryError.getStackTrace().length > 0)
                Log.error("Out of memory stack-trace: "+outOfMemoryError.getStackTrace()[0].toString(), outOfMemoryError);
            System.exit(-1);
        } catch (IOException ioError) {
            // check for some specific errors
            if (ioError.getMessage().equals("No space left on device")) {
                Log.error("[DISK] There is no space left on the device!");
            } else {
                Log.error("Error while executing " + tool.getClass(), ioError);
            }
            System.exit(-1);
        } catch (Exception e) {
            Log.error("","\n");
            if(e.getMessage() != null)
                Log.error(e.getMessage(), e);   // always provide stacktrace
            else
                Log.error(e.getMessage(), e);
            Log.error("","");
            Log.debug("\n\n");
            Log.debug("Error while executing " + tool.getClass() + " : " + e.getMessage(), e);
            System.exit(-1);
        } finally {
            // shutdown the executor
            Execute.shutdown();
        }
    }

    /**
     * Create the default flux command options
     * @return
     */
    private static JSAP createBaseOptions() {
        /*
       Create command line parser and
       add default Flux parameter
        */
        JSAP jsap = new JSAP();
        try {
            jsap.registerParameter(JSAPParameters.flaggedParameter("tool", 't').defaultValue(System.getProperty("flux.tool")).help("Select a tool").get());
            jsap.registerParameter(JSAPParameters.switchParameter("help", 'h').help("Show help").get());
            jsap.registerParameter(JSAPParameters.switchParameter("list-tools").help("List available tools").get());
            jsap.registerParameter(JSAPParameters.flaggedParameter("threads").defaultValue("2").type(Integer.class).help("Maximum number of threads to use. Default 2").get());
            jsap.registerParameter(JSAPParameters.flaggedParameter("log").defaultValue("INFO").help("Log level (NONE|INFO|ERROR|DEBUG)").valueName("level").get());
            jsap.registerParameter(JSAPParameters.switchParameter("force").help("Disable interactivity. No questions will be asked").get());
            jsap.registerParameter(JSAPParameters.switchParameter("version", 'v').help("Show version information").get());
        } catch (JSAPException e) {
            Log.error("Unable to create parameters : " + e.getMessage(), e);
            System.exit(-1);
        }
        return jsap;
    }

    /**
     * Sets a <code>SecurityManager</code> when required by property
     * variables.
     * @return a security manager or <code>null</code> if none is needed
     */
    protected static SecurityManager getSecurityManager() {
    	
    	final Properties props= System.getProperties();
    	Iterator iter= props.keySet().iterator();
    	boolean createSM= false;
    	for(;iter.hasNext()&& !createSM;) 
    		if (iter.next().toString().startsWith("flux.io.deny"))
    			createSM= true;
    	if (!createSM)
    		return null;
    	
		SecurityManager sm= new SecurityManager() {
			File systemTemp= new File(System.getProperty("java.io.tmpdir"));
			@Override
			public void checkPermission(Permission perm) {
				checkPermission(perm, null);
			};
			@Override
			public void checkPermission(Permission perm, Object context) {
				if (perm instanceof FilePermission) {
					if (props.containsKey("flux.io.deny.tmpdir")
								&& (perm.getActions().contains("write"))) {
							File f= new File(perm.getName());
							if (!(f.exists()|| f.isDirectory()))
								f= f.getParentFile();
							if (f.equals(systemTemp))
								throw new AccessControlException("access denied", perm);
					}
				}
			}
		};

		return sm;
	}

	public static void printUsage(Tool tool, JSAP jsap, List<Tool> allTools, String errorMessage, boolean userRequestsHelp) {
        if(errorMessage != null){
            System.err.println(errorMessage);
            System.err.println("");
            System.exit(1);
        }

        if(userRequestsHelp){
            System.err.println("-------Documentation & Issue Tracker-------");
            System.err.println("Flux Wiki (Docs): http://sammeth.net/confluence");
            System.err.println("Flux JIRA (Bugs): http://sammeth.net/jira");
            System.err.println("");
            System.err.println("Please feel free to create an account in the public");
            System.err.println("JIRA and reports any bugs or feature requests.");
            System.err.println("-------------------------------------------");
            System.err.println("");

            if(tool != null){
                System.err.println("Current tool: " + tool.getName());
                System.err.println("");
                String description = tool.getLongDescription();
                if(description == null) description = tool.getDescription();
                if(description != null){
                    System.err.println(description);
                    System.err.println("");
                }
                // custom jsap for the tool
                JSAP toolJSAP = new JSAP();
                List<Parameter> parameter = tool.getParameter();
                if(parameter != null){
                    try{
                        for (Parameter p : parameter) {
                            toolJSAP.registerParameter(p);
                        }
                    } catch (Exception e) {
                        Log.error("Parameter error : " + e.getMessage(), e);
                        System.exit(-1);
                    }
                }

                System.err.println("Tool specific options:\n");
                System.err.println(toolJSAP.getHelp());
                System.err.println("");
            }

            JSAP baseOptions = createBaseOptions();
            System.err.println("The Flux library comes with a set of tools.\n" +
                    "You can switch tools with the -t option. The general options\n" +
                    "change the behaviour of all the packaged tools.");
            System.err.println("");
            System.err.println("General Flux Options: \n");
            System.err.println(baseOptions.getHelp());
            System.err.println("");

            printTools(allTools);
            System.err.flush();
        }
        System.exit(-1);
    }

    /**
     * Print available flux tools
     *
     * @param tools the tools available flux tools
     */
    private static void printTools(List<Tool> tools) {
        System.err.println("The Flux library consists of a set of tools bundled with the package.");
        if(System.getProperty("flux.tool", null) != null){
            System.err.println("The current bundle uses '" + System.getProperty("flux.tool") + "' as the default tool.");
        }
        System.err.println("You can switch tools with the -t option and get help for a specific\n" +
                "tool with -t <toolname> --help. This will print the usage and description of the specified tool");

        System.err.println("");
        System.err.println("List of available tools");
        for (Tool fluxTool : tools) {
            System.err.println("\t" + fluxTool.getName() + " - " + fluxTool.getDescription());
        }
        System.err.println();
    }

    /**
     * Search for flux tools. Currently the search is limited to classes in the barna package.
     *
     * @return tools all detected flux tools
     */
    public static List<Tool> findTools() {

        // scan the classpath to find tools
        List<Tool> tools = new ArrayList<Tool>();
        List<URL> fbiUrls = null;
        try {
        	fbiUrls= new ArrayList<URL>(ClasspathHelper.forPackage("barna"));
        } catch (Throwable t) {
        	System.err.println(t);
        	System.currentTimeMillis();
        }

        ConfigurationBuilder config = new ConfigurationBuilder();
        config.useParallelExecutor();
        config.setScanners(new SubTypesScanner());
        config.setUrls(fbiUrls);

        Reflections reflections = new Reflections(config);


        Set<Class<? extends Tool>> toolClasses = reflections.getSubTypesOf(Tool.class);
        for (Class<? extends Tool> toolClass : toolClasses) {
            try {
                // BARNA-304 -- added check for interface extensions or
                // abstract class implementation of Tool
                if(!toolClass.isInterface() && !Modifier.isAbstract(toolClass.getModifiers())){
                    Tool fluxTool = toolClass.newInstance();
                    tools.add(fluxTool);
                }
            } catch (Exception e) {
                Log.error("Error while creating tool instance for " + toolClass.getName(), e);
                Log.error("Make sure the class exists and has a default constructor.");
            }
        }
        return tools;
    }

    /**
     * Get the name of the tool that is used
     *
     * @return toolName the name of the tool
     */
    public String getToolName() {
        return toolName;
    }

    /**
     * Set the name of the tool
     *
     * @param toolName toolname
     */
    public void setToolName(String toolName) {
        this.toolName = toolName;
    }


    /**
     * Set the log level. The Sting must be one of {@code NONE, INFO, ERROR, DEBUG}
     *
     * @param level the level
     */
    public void setLogLevel(String level) {
        Log.setLogLevel(level);
    }

    /**
     * Get the current log level
     *
     * @return level the log level
     */
    public String getLogLevel() {
        return Log.getLogLevel().toString();
    }

    /**
     * Disable interactivity
     *
     * @param detached detache
     */
    public void setDetached(boolean detached) {
        Log.setInteractive(!detached);
    }

    /**
     * Returns true if interactivity is disabled
     *
     * @return detached true if interactivity is disabled
     */
    public boolean isDetached() {
        return !Log.isInteractive();
    }

    /**
     * Returns the number of available background threads
     *
     * @return threads the number of background threads
     */
    public int getThreads() {
        return threads;
    }

    /**
     * Set the number of background threads
     *
     * @param threads number of background threads
     */
    public void setThreads(final int threads) {
        this.threads = threads;
    }

    /**
     * Read properties like version and build revision from jar file
     *
     * @return valid returns true if properties are valid and everything is fine
     */
    private static void checkJavaVersion() {
        try {
            float v = FLUX_JAVA_VERSION;
            String ver = System.getProperty("java.version");
            int p = ver.indexOf('.', 0);
            p = ver.indexOf('.', p + 1);
            float v2 = Float.parseFloat(ver.substring(0, p));
            if (v2 < v) {
                Log.error("Wrong java version, I need " + v + " but I found " + v2 + ".");

            }
        } catch (Exception e) {
            ; // :)
        }
    }

    /**
     * Cover build and version information
     */
    public static class FluxVersionInfo{
        private String appName;
        private String libVersion;
        private String appVersion;
        private String buildDate;
        private String buildVersion;
        private String buildBranch;

        /**
         * Createa a new version info from a properties file. File path is relative to the classpath, Flux
         * is used to find the resource.
         *
         * @param fileName the file name (class path reference)
         */
        public FluxVersionInfo(String fileName){
            InputStream buildProperties = getClass().getResourceAsStream(fileName);
            if(buildProperties != null){
                Properties properties = new Properties();
                try {
                    properties.load(buildProperties);
                    libVersion = properties.getProperty("barna.version", "Unknown");
                    appVersion = properties.getProperty("flux.appversion", "Unknown");
                    appName = properties.getProperty("flux.appname", "Flux");
                    buildDate = properties.getProperty("build.date", "Unknown");
                    buildBranch = properties.getProperty("build.branch", "Unknown");
                    buildVersion = properties.getProperty("build.version", "Unknown");
                    buildVersion = properties.getProperty("build.version", "Unknown");
                } catch (IOException ignore) {
                    // ignore
                }
            }
        }

        /**
         * Get the flux library version
         * @return libVersion the library version
         */
        public String getLibVersion() {
            return libVersion;
        }

        /**
         * Get the application version
         * @return appVersion the application version
         */
        public String getAppVersion() {
            return appVersion;
        }

        /**
         * Get the build date
         * @return buildDate the build date
         */
        public String getBuildDate() {
            return buildDate;
        }

        /**
         * Get the build version
         * @return buildVersion the build version
         */
        public String getBuildVersion() {
            return buildVersion;
        }

        /**
         * Get the branch name
         * @return buildBranch the branch name
         */
        public String getBuildBranch() {
            return buildBranch;
        }

        /**
         * Get the application name
         * @return name the name
         */
        public String getAppName() {
            return appName;
        }

        /**
         * Long multi-line string representation of all the version information
         * @return info the build info
         */
        public String toString(){
            StringBuilder builder = new StringBuilder();
            builder.append(appName).append('\n');
            builder.append("Version ").append(appVersion).append('\n');
            builder.append("Barna Library ").append(libVersion).append('\n');
            builder.append("-----------------------------------------------\n");
            builder.append("Build Date ").append(buildDate).append('\n');
            builder.append("Build Version ").append(buildVersion).append('\n');
            builder.append("Build Branch ").append(buildBranch).append('\n');
            return builder.toString();
        }
        /**
         * Shor version info
         * @return info short version info
         */
        public String toShortString(){
            StringBuilder builder = new StringBuilder();
            builder.append(appName).append(" ").append("v").append(appVersion).append(" (Flux Library: ").append(libVersion).append(")\n");
            return builder.toString();
        }
    }
}
