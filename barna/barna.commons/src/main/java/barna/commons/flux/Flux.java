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

package barna.commons.flux;

import barna.commons.Execute;
import barna.commons.Log;
import barna.commons.options.HelpPrinter;
import org.cyclopsgroup.caff.ref.AccessFailureException;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;
import org.cyclopsgroup.jcli.spi.ParsingContext;
import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.io.File;
import java.io.FilePermission;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.security.AccessControlException;
import java.security.Permission;
import java.util.*;


/**
 * Flux Simulator starter class. Contains the main method and parses command line arguments. During startup,
 * this checks for any {@link FluxTool} implementations and adds them to the list of available tools.
 */
@Cli(name = "flux", restrict = false)
public class Flux {
    /**
     * Current Flux Simulator version
     */
    public static String FLUX_VERSION = "";

    /**
     * Current Flux Simulator revision
     */
    public static String FLUX_REVISION = "";

    /**
     * The tool to be executed
     */
    private String toolName;

    /**
     * Show help message
     */
    private boolean help;

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
        /*
        Read properties
         */
        if (!readProperties()) {
            System.exit(-1);
        }

        // register tools
        List<FluxTool> tools = findTools();

        // prepare the simulator
        Flux simulator = new Flux();
        ArgumentProcessor fluxArguments = ArgumentProcessor.newInstance(Flux.class);
        try {
            fluxArguments.process(args, simulator);
        } catch (AccessFailureException ae) {
            Log.error("Error while processing arguments !");
            if (ae.getCause() instanceof InvocationTargetException) {
                Log.error(((InvocationTargetException) (ae.getCause())).getTargetException().getMessage());
            } else {
                Log.error(ae.getCause().getMessage());
            }
            System.exit(-1);
        } catch (Exception e) {
            Log.error("Error while processing arguments !");
            Log.error(e.getMessage());
            System.exit(-1);
        }


        // prepare tools
        List<org.cyclopsgroup.jcli.spi.Cli> toolClis = new ArrayList<org.cyclopsgroup.jcli.spi.Cli>();
        for (FluxTool fluxTool : tools) {
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
            ParsingContext context = toolArguments.createParsingContext();
            toolClis.add(context.cli());
        }


        // start
        if (FLUX_VERSION.length() > 0 && FLUX_REVISION.length() > 0) {
            Log.info("I am the Flux Toolbox (v" + FLUX_VERSION + " build " + FLUX_REVISION + "), nice to meet you!\n");
        } else {
            Log.info("I am the Flux Toolbox ( Devel Mode ), nice to meet you!\n");
        }

        // find the tool to start
        FluxTool tool = null;
        if (simulator.getToolName() != null) {
            int i = 0;
            for (org.cyclopsgroup.jcli.spi.Cli cli : toolClis) {
                if (cli.getName().equals(simulator.getToolName())) {
                    tool = tools.get(i);
                    break;
                }
                i++;
            }
        }

        // delegate help prints to std err
        PrintWriter out = new PrintWriter(System.err);
        HelpPrinter printer = new HelpPrinter(out);

        if (simulator.isHelp()) {
            // show help message
            if (tool == null) {
                printFluxHelp(tools, fluxArguments, out, printer);
            } else {
                // show tool help
                // create the argument parser
                ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
                //toolArguments.printHelp(out);
                printer.print(toolArguments);
                out.flush();
            }
            // exit after printing help
            System.exit(-1);
        }


        if (tool != null) {
            // execute the tool
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
            toolArguments.process(args, tool);
            if (!tool.validateParameters(printer, toolArguments)) {
                out.flush();
                System.exit(-1);
            }

            try {
                // configure the executor
                Execute.initialize(simulator.getThreads());

                tool.call();
            } catch (IOException ioError) {
                // check for some specific errors
                if (ioError.getMessage().equals("No space left on device")) {
                    Log.error("[DISK] There is no space left on the device!");
                } else {
                    Log.error("Error while executing " + tool.getClass(), ioError);
                }
                System.exit(-1);
            } catch (Exception e) {
                Log.error(e.getMessage());
                Log.debug("\n\n");
                Log.debug("Error while executing " + tool.getClass() + " : " + e.getMessage(), e);
                System.exit(-1);
            } finally {
                // shutdown the executor
                Execute.shutdown();
            }
        } else {
            if (simulator.getToolName() == null || simulator.getToolName().isEmpty()) {
                Log.error("");
                Log.error("No tool specified!");
                Log.error("\n");
            } else {
                Log.error("");
                Log.error("Unable to find tool : " + simulator.getToolName());
                Log.error("\n");
            }

            // warn and show help
            printFluxHelp(tools, fluxArguments, out, printer);
            System.exit(-1);
        }
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

	private static void printFluxHelp(List<FluxTool> tools, ArgumentProcessor fluxArguments, PrintWriter out, HelpPrinter printer) {
        // show general help
        //fluxArguments.printHelp(out);
        printer.print(fluxArguments);
        out.println("\tTools available:");
        for (FluxTool fluxTool : tools) {
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
            ParsingContext context = toolArguments.createParsingContext();
            out.println("\t\t" + context.cli().getName() + " - " + context.cli().getDescription());
        }
        out.println();
        out.println("To get help for a specific tool try -t <tool> --help");
        out.println();
        out.flush();
    }

    /**
     * Search for flux tools. Currently the search is limited to classes in the barna package.
     *
     * @return tools all detected flux tools
     */
    public static List<FluxTool> findTools() {

        // scan the classpath to find tools
        List<FluxTool> tools = new ArrayList<FluxTool>();
        List<URL> fbiUrls = null;
        try {
        	fbiUrls= new ArrayList<URL>(ClasspathHelper.getUrlsForPackagePrefix("barna"));
        } catch (Throwable t) {
        	System.err.println(t);
        	System.currentTimeMillis();
        }

        ConfigurationBuilder config = new ConfigurationBuilder();
        config.useParallelExecutor();
        config.setScanners(new SubTypesScanner());
        config.setUrls(fbiUrls);

        Reflections reflections = new Reflections(config);


        Set<Class<? extends FluxTool>> toolClasses = reflections.getSubTypesOf(FluxTool.class);
        for (Class<? extends FluxTool> toolClass : toolClasses) {
            try {
                FluxTool fluxTool = toolClass.newInstance();
                tools.add(fluxTool);
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
    @Option(name = "t", longName = "tool", description = "select the tool", displayName = "tool", required = true)
    public void setToolName(String toolName) {
        this.toolName = toolName;
    }

    /**
     * Should the help be displayed
     *
     * @return help true if help message should be displayed
     */
    public boolean isHelp() {
        return help;
    }

    /**
     * Show the help message
     *
     * @param help the help message
     */
    @Option(name = "h", longName = "help", description = "show help message")
    public void setHelp(boolean help) {
        this.help = help;
    }

    /**
     * Set the log level. The Sting must be one of {@code NONE, INFO, ERROR, DEBUG}
     *
     * @param level the level
     */
    @Option(name = "", longName = "log", description = "Log level (NONE|INFO|ERROR|DEBUG)", defaultValue = "INFO", required = false)
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
    @Option(name = "f", longName = "force", description = "Disable interactivity. No questions will be asked", required = false)
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
    @Option(name = "x", longName = "threads", description = "Number of threads, default is 2", required = false, defaultValue = "2")
    public void setThreads(final int threads) {
        this.threads = threads;
    }

    /**
     * Read properties like version and build revision from jar file
     *
     * @return valid returns true if properties are valid and everything is fine
     */
    private static boolean readProperties() {
        /*
        Find the manifest file and extract version revision adn jdk information
         */
        URL location = Flux.class.getResource("/flux.properties");
        String buildVersion = "";
        String buildRevision = "";
        String buildJDK = "";
        if (location != null) {
            try {
                Properties pp = new Properties();
                pp.load(location.openStream());
                String v = pp.getProperty("Build-Version");
                String r = pp.getProperty("buildNumber");
                String j = pp.getProperty("Flux-JDK");
                if (v != null) {
                    buildVersion = v.toString();
                }
                if (r != null) {
                    buildRevision = r.toString();
                }
                if (j != null) {
                    buildJDK = j.toString();
                }
            } catch (IOException e) {
            }
        }


        if (!buildVersion.isEmpty()) {
            FLUX_VERSION = buildVersion;
        }

        if (!buildRevision.isEmpty()) {
            FLUX_REVISION = buildRevision;
        }

        if (!buildJDK.isEmpty()) {
            try {
                float v = Float.parseFloat(buildJDK);
                String ver = System.getProperty("java.version");
                int p = ver.indexOf('.', 0);
                p = ver.indexOf('.', p + 1);
                float v2 = Float.parseFloat(ver.substring(0, p));
                if (v2 < v) {
                    Log.error("Wrong java version, I need " + v + " but I found " + v2 + ".");
                    return false;
                }
            } catch (Exception e) {
                ; // :)
            }

        }
        return true;

    }
}
