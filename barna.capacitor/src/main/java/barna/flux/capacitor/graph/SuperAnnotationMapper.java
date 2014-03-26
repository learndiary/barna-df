package barna.flux.capacitor.graph;

import barna.model.SuperLocus;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/7/13
 * Time: 6:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class SuperAnnotationMapper extends AnnotationMapper {

    AnnotationMapper[] annos= null;

    public SuperAnnotationMapper(SuperLocus sl) {
        //super(null, null, false);
        super(sl, false, false, false, null);
    }

}
