package fbi.commons.io;

/**
 * Use this to access the IOHandler implementations
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class IOHandlerFactory {

    /**
     * No instances allowed
     */
    private IOHandlerFactory(){}

    /**
     * THe default handler
     */
    private static IOHandler defaultHandler;

    /**
     * Get the default handler
     *
     * @return ioHandler the default handler
     */
    public static IOHandler getDefaultHandler(){
        if(defaultHandler == null){
            defaultHandler = new SimpleIOHandler();
        }
        return defaultHandler;
    }
}
