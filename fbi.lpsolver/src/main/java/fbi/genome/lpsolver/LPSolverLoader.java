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

package fbi.genome.lpsolver;

import lpsolve.LpSolve;
import lpsolve.VersionInfo;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Check the operating system and load the native libraries
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class LPSolverLoader {
    /**
     * Library name on linux
     */
    private static String FILENAME_LX="liblpsolve55.so";
    /**
     * Library name on linux for JNI binding
     */
    private static String FILENAMEJ_LX="liblpsolve55j.so";
    /**
     * Library name on OSX
     */
    private static String FILENAME_OSX="liblpsolve55.dylib";
    /**
     * Library name on OSX for JNI binding
     */
    private static String FILENAMEJ_OSX="liblpsolve55j.jnilib";
    /**
     * Library name on Win
     */
    private static String FILENAME_WIN="lpsolve55.dll";
    /**
     * Library name on Win for JNI binding
     */
    private static String FILENAMEJ_WIN="lpsolve55j.dll";

    /**
     * Load the library
     *
     * @return success true if loaded
     */
    public static boolean load(){
        String bits = OSChecker.is64bit() ? "64":"32";

        String name = null;
        String jname = null;
        String dir = null;
        if(OSChecker.isWindows()){
            dir = "win";
            name = FILENAME_WIN;
            jname = FILENAMEJ_WIN;
        }else if(OSChecker.isMacOSX()){
            dir = "osx";
            name = FILENAME_OSX;
            jname = FILENAMEJ_OSX;
        }else if(OSChecker.isLinux()){
            dir = "linux";
            name = FILENAME_LX;
            jname = FILENAMEJ_LX;
        }
        if(dir == null){
            throw new RuntimeException("LPSolve library not found for your operating system!");
        }

        // check if the files exist, otherwise copy the lib from jar to file
        File tmpDir = new File(System.getProperty("java.io.tmpdir"));
        File libFile = new File(tmpDir, name);
        File libjFile = new File(tmpDir, jname);

        //if(!libFile.exists()){
            // copy
            try {
                URL dylib = LPSolverLoader.class.getResource("/" + dir + "/" + bits + "/" + name);
                if(!OSChecker.isMac() && !OSChecker.isLinux()){ // mac os is a single jnilib file
                    write2File(dylib, libFile);
                }
            } catch (IOException e) {
                throw new RuntimeException("Unable to copy library to filesystem! " + e.getMessage(), e);
            }
        //}
        //if(!libjFile.exists()){
            // copy
            try {
                write2File(LPSolverLoader.class.getResource("/"+dir+"/"+bits+"/"+jname), libjFile);
            } catch (IOException e) {
                throw new RuntimeException("Unable to copy library to filesystem! " + e.getMessage(), e);
            }
        //}

        // files should exist now, load
        /*
        CAP-7 Make sure we only load the shared library on windows
        as we have single JNI library for OSX and LInux
         */
        if(libFile.exists() && !OSChecker.isMac() && !OSChecker.isLinux()){
            System.load(libFile.getAbsolutePath());
        }

        if(libjFile.exists()){
            System.load(libjFile.getAbsolutePath());
        }

        // check
        VersionInfo versionInfo = LpSolve.lpSolveVersion();
        return versionInfo != null;
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
            throw new RuntimeException("unable to write URL to Directory " + file.getAbsolutePath());
        }
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
    }

}
