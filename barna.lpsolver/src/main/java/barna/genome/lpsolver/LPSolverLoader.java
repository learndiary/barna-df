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

package barna.genome.lpsolver;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Supports loading the LPSolver native libraries.
 * The loader supports two environment variables:
 * <ul>
 *     <li>LPSOLVER_JNI - for the JNI library</li>
 *     <li>LPSOLVER_LIB - for the shared object library</li>
 * </ul>,
 * You can specify the environment variables to load from a given location
 * otherwise the bundled libraries are copied to the users home directory.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class LPSolverLoader {
    /**
     * Indicate load status for this VM
     */
    private static boolean IS_LOADED = false;
    /**
     * If load failed once, we set this to true
     */
    private static boolean TRIED_LOADING = false;

    /**
     * Library name on linux
     */
    private static final String FILENAME_LX="liblpsolve55.so";
    /**
     * Library name on linux for JNI binding
     */
    private static final String FILENAMEJ_LX="liblpsolve55j.so";
    /**
     * Library name on OSX for JNI binding
     */
    private static final String FILENAMEJ_OSX="liblpsolve55j.jnilib";
    /**
     * Library name on Win
     */
    private static final String FILENAME_WIN="lpsolve55.dll";
    /**
     * Library name on Win for JNI binding
     */
    private static final String FILENAMEJ_WIN="lpsolve55j.dll";
    /**
     * Environment variable checked for JNI library file
     */
    public static final String ENV_JNI = "LPSOLVER_JNI";
    /**
     * Environment variable checked for library file
     */
    public static final String ENV_LIB = "LPSOLVER_LIB";

    /**
     * Load the library. This checks for environment variables
     * and eventually copies the bundled libraries
     *
     * @throws IOException in case of any copy errors
     * @throws UnsatisfiedLinkError in case the library loading failed
     */
    public static void load() throws IOException, UnsatisfiedLinkError {
        if(TRIED_LOADING) return;
        TRIED_LOADING = true;
        /*
        CAP-8 try to write to user home
         */
        File targetDirectory = new File(new File(System.getProperty("user.home"), ".lpsolver"), OSChecker.is32bit() ? "lpsolve_32": "lpsolve_64");
        if(!targetDirectory.exists()){
            if(!targetDirectory.mkdirs()){
                // unable to write to user home ... for whatever reason
                // switch to global tmp
                targetDirectory = new File(System.getProperty("java.io.tmpdir"), OSChecker.is32bit() ? "lpsolve_32": "lpsolve_64");
                targetDirectory.mkdirs();
            }
        }

        // prepare the libraries
        File[] libFiles = prepareLibraries(targetDirectory);
        load(libFiles[0], libFiles[1]);
    }

    /**
     * Checks for environment variables and eventually copies bundled libs
     * to target directory if they do not exist yet
     *
     * @param targetDirectory the target directory
     * @return files returns an array, first element is the jniLibrary file, second is the shared library which can also be null
     */
    private static File[] prepareLibraries(File targetDirectory) throws IOException {
        File[] files = new File[2];
        // set up the classpath directory
        String dir = getClassPathDirectory();

        // first check for the JNI library
        if(System.getenv(ENV_JNI) != null){
            files[0] = new File(System.getenv(ENV_JNI));
        }else{
            // get a filename
            String jniFile = getJNILibraryName();
            files[0] = new File(targetDirectory, jniFile);
            if(!files[0].exists()){
                URL urlToLib = LPSolverLoader.class.getResource(dir + jniFile);
                write2File(urlToLib, files[0]);
            }
        }


        // second check for the shared library
        if(System.getenv(ENV_LIB) != null){
            files[1] = new File(System.getenv(ENV_LIB));
        }else{
            // get a filename
            String libFile = getSharedLibraryName();
            if(libFile != null){
                files[1] = new File(targetDirectory, libFile);
                if(!files[1].exists()){
                    URL urlToLib = LPSolverLoader.class.getResource(dir + libFile);
                    write2File(urlToLib, files[1]);
                }
            }
        }

        return files;
    }

    /**
     * Created the path to the libraries folder in the classpath to install
     * bundled libraries. NOTE that you have to append the filename !
     *
     * @return classpathDir the path to the dir - you still have to append the filename
     */
    private static String getClassPathDirectory() {
        String dir = "/";
        if(OSChecker.isWindows()){
            dir += "win/";
        }else if(OSChecker.isMacOSX()){
            dir += "osx/";
        }else{
            dir += "linux/";
        }

        if(OSChecker.is64bit()){
            dir += "64/";
        }else{
            dir += "32/";
        }
        return dir;
    }

    /**
     * Load the library files manually. The shared library can be null
     *
     * @param jniLibrary the native interface
     * @param sharedLibrary the shared library
     * @throws Exception in case the library could not be loaded successfully
     */
    public static void load(File jniLibrary, File sharedLibrary) throws UnsatisfiedLinkError{
        if(jniLibrary == null ) throw new NullPointerException("You have to specify a JNI library file");
        if(!jniLibrary.exists()) throw new IllegalArgumentException("The specified JNI library " + jniLibrary + " does not exist");
        if(sharedLibrary != null && !sharedLibrary.exists()) throw new IllegalArgumentException("The specified shared library " + sharedLibrary + " does not exist");
        // avoid loading twice
        if(TRIED_LOADING && IS_LOADED) return;

        try{
            TRIED_LOADING = true;
            if(sharedLibrary != null){
                System.load(sharedLibrary.getAbsolutePath());
            }
            System.load(jniLibrary.getAbsolutePath());
            IS_LOADED = true;
        }catch(UnsatisfiedLinkError error){
            getLogger().log(Level.SEVERE, "Error while loading the lpsolver native libraries: " + error.getMessage(), error);
            throw error;
        }
    }

    /**
     * Helper to write URL content to file
     *
     * @param url the source URL
     * @param file the target file
     * @throws IOException in case of any errors
     */
    private static void write2File(URL url, File file) throws IOException {
        if(file.isDirectory()){
            throw new RuntimeException("Unable to write URL to Directory " + file.getAbsolutePath());
        }
        try{
            getLogger().info("Copy library to " + file.getAbsolutePath());
            if(!file.exists()){
                file.createNewFile();
            }
            FileOutputStream out = null;
            InputStream in = null;
            try {
                in = url.openStream();
                out = new FileOutputStream(file);
                byte[] buff = new byte[4096];
                int len = 0;
                while (-1 != (len = in.read(buff))) {
                    out.write(buff, 0, len);
                }
            } finally{
                if(out != null)out.close();
                if(in != null) in.close();
            }
        }catch(IOException error){
            getLogger().log(Level.SEVERE, "Unable to copy bundled library to " + file.getAbsolutePath(), error);
            throw error;
        }
    }

    private static Logger getLogger() {
        return Logger.getLogger(LPSolverLoader.class.getName());
    }

    /**
     * Get the name of the JNI library file depending on the current OS
     *
     * @return name the file name
     * @throws RuntimeException in case the OS is not supported
     */
    private static String getJNILibraryName(){
        if(OSChecker.isWindows()){
            return FILENAMEJ_WIN;
        }else if(OSChecker.isMacOSX()){
            return FILENAMEJ_OSX;
        }else if(OSChecker.isLinux()){
            return FILENAMEJ_LX;
        }
        throw new RuntimeException("The current Operating system is not supported!");
    }

    /**
     * Get the name of the shared library file depending on the current OS
     *
     * @return name the file name or null if no shared library exists for the current OS
     */
    private static String getSharedLibraryName(){
        if(OSChecker.isWindows()){
            return FILENAME_WIN;
        }else if(OSChecker.isLinux()){
            return FILENAME_LX;
        }
        return null;
    }

}
