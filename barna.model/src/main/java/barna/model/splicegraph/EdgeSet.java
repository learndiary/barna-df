package barna.model.splicegraph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Administrates a set of edges.
 */
public class EdgeSet {

    /**
     * Set to handle unique instances for each edge group.
     */
    public HashMap<EdgeSet, EdgeSet> singletonMgr;   // do not make static, will accumulate across loci!

    /**
     * Returns the set of edges managed by <code>this</code> instance.
     * @return a hashset with the edges
     */
    public HashMap<AbstractEdge, Byte> getEset() {
        return eset;
    }

    HashMap<AbstractEdge, Byte> eset; // TODO replace by constant growth hash, possibly on disk

    // for internal use only
    private EdgeSet() {
    }

    public EdgeSet(AbstractEdge e, byte dir, int size) {
        singletonMgr= new HashMap<EdgeSet, EdgeSet>(100, 1f);
        eset= new HashMap<AbstractEdge, Byte>(size, 1f);
        add(e, dir);
    }

    /**
     * Adds non-redundantly an element to the edge set.
     * @param e the edge to be added
     * @return <code>this</code> regardless whether the edge was added to the current set or not,
     * or another instance if the extended set already exists in the collection of singletons
     */
    public EdgeSet add(AbstractEdge e, Byte dir){

        if (eset.containsKey(e)&& eset.get(e)== dir)
            return this;

        eset.put(e, dir);
        if (singletonMgr.containsKey(e)) {
            eset.remove(e); // in case the object is stored as singleton
            return singletonMgr.get(e);
        }

        // else
        singletonMgr.put(this, this);
        return this;
    }

    @Override
    public int hashCode() {

        StringBuilder sb= new StringBuilder();
        for (AbstractEdge e:eset.keySet()) {    // disregard directionalities
             sb.append(e.toString());
        }

        return sb.toString().hashCode();
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof EdgeSet))
            return false;
        EdgeSet otherEset= (EdgeSet) obj;

        if (eset.size()!= otherEset.eset.size())
            return false;

        // falls back to edge and ss comparison, incl chr
        Iterator<AbstractEdge> iter= eset.keySet().iterator();
        while(iter.hasNext()) {
            AbstractEdge k= iter.next();
            if ((!otherEset.eset.containsKey(k))
                || (otherEset.eset.get(k)!= eset.get(k)))
                return false;
        }

        return true;
    }


    /**
     * Obtains the current size of the set.
     * @return size of the underlying set of edges.
     */
    public int size() {
        return eset.size();
    }

}
