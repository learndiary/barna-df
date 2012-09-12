package barna.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 5/9/12
 * Time: 12:22 PM
 */
public interface MSIteratorDisk<T> extends MSIterator<T>{

    /**
     * Obtains the sorted file the <code>BEDiteratorDisk</code> instance
     * is based on. <b>Note</b>: intrinsically invokes <code>init()</code>
     * if no explicit call has incurred yet.
     * @return the sorted file on which <code>this</code> iterator is
     * based on
     * @throws java.util.concurrent.ExecutionException
     * @throws IOException
     * @see #init()
     */
    public File getTmpFile() throws ExecutionException, IOException;

    /**
     * Produces a sorted file for subsequent iteration from the
     * correspondingly provided input (i.e., streams or files,
     * sorted or unsorted). Possible copying and sorting processes
     * are executed in parallel, It is recommended to call <code>init()</code>
     * directly after the constructor. Otherwise, an intrinsic
     * call will ensure correct initialization by the first time
     * <code>hasNext()</code> or <code>next()</code> are invoked.
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     * @see #getTmpFile()
     */
    public void init() throws FileNotFoundException, IOException;

    /**
     * Determines the number of elements left in the
     * iterator, starting to count at the current
     * position.<br>
     * <b>Warning</b>: time complexity of the method
     * scales linearly with the elements that are
     * left and have to be read from the file.
     * @return the number of elements from the current
     * reading position until the end
     */
    public int countRemainingElements();

}
