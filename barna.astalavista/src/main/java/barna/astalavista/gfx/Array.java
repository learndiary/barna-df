package barna.astalavista.gfx;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/6/12
 * Time: 11:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class Array {
    Object[] array= null;

    public Array(Object[] a) {
        array= a;
    }

    @Override
    public int hashCode() {
        if (array== null|| array.length< 1)
            return 0;
        return array[0].hashCode();	// assuming object (not value!) identity of elements
    }

    @Override
    public boolean equals(Object obj) {
        if (obj== null|| !(obj instanceof Array))
            return false;
        Array a= (Array) obj;
        if (array== null^a.array== null)
            return false;
        if (array== null&&a.array== null)
            return true;
        if (array.length!= a.array.length)
            return false;
        for (int i = 0; i < array.length; i++)
            if (array[i]!= a.array[i])
                return false;
        return true;
    }
}
