package barna.io;

import barna.commons.parameters.FileNameParser;

import java.io.File;

/**
 * Implement a file name parser that works relative to a given parent directory
 */
public class RelativePathParser implements FileNameParser {
    /**
     * the Parent directory
     */
    File parentDir = null;

    @Override
    public File parse(String string) {
        return FileHelper.fromRelative(string, parentDir);
    }

    /**
     * Get the current parent directory
     *
     * @return parent the current parent
     */
    public File getParentDir() {
        return parentDir;
    }

    /**
     * Set the current parent directory
     *
     * @param parentDir the parent directory
     */
    public void setParentDir(File parentDir) {
        this.parentDir = parentDir;
    }
}
