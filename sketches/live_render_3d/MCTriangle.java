package rui.marchingCubes;

import toxi.geom.Vec3D;

/**
 * simple container for triangle vertices and normal
 * @author ruimadeira
 *
 */

public class MCTriangle {
	public Vec3D a, b, c, normal;
	
	MCTriangle(){
		
	}
	MCTriangle(Vec3D _a, Vec3D _b, Vec3D _c){
		a = new Vec3D(_a);
		b = new Vec3D(_b);
		c = new Vec3D(_c);
		normal = new Vec3D();
	}
	MCTriangle(Vec3D _a, Vec3D _b, Vec3D _c, Vec3D _norm){
		a = new Vec3D(_a);
		b = new Vec3D(_b);
		c = new Vec3D(_c);
		normal = new Vec3D(_norm);
	}
}
