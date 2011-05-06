package fbi.commons.tools;

/**
 * Implementations deal with object and support intercepting the object at specific points. For example,
 * the {@link StreamSorter} implementation can implement this interface to allow the user to intercept
 * the sorted lines without dealing with a separate thread and piped streams.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface Interceptable<T> {
    void addInterceptor(Interceptor<T> interceptor);

    /**
     * Interceptors catch object from the Interceptable. They can return a changed version of the
     * line, but it is up to the Intercaptable to accept the change.
     */
    public static interface Interceptor<T> {
        /**
         * Intercepts an object. Implementations can return a changed version but it is up to the
         * interceptable to accept the change
         *
         * @param object the intercepted object
         * @return object the mutated object
         */
        T intercept(T object);
    }
}
