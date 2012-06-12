package barna.model.splicegraph;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 6/12/12
 * Time: 10:25 AM
 */
public class AllIntronicEdge extends SimpleEdge {

    public AllIntronicEdge(Node newTail, Node newHead) {
        super(newTail, newHead);
    }

    @Override
    public boolean isExonic() {
        return false;
    }

    @Override
    public boolean isIntronic() {
        return false;
    }

    @Override
    public boolean isAllIntronic() {
        return true;
    }
}
