package rui.marchingCubes;

import processing.core.PApplet;

import java.lang.Math;
import java.util.ArrayList;
import toxi.geom.Vec3D;

/**
 * 
 * 
 * simple class implementing the Marching Cubes algorithm to create 3d volumetric meshes.
 * based on the code and explanations by Paul Bourke that can be found here:
 * http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
 * 
 * its dependent on Processing's PApplet and Karsten Schmidt's Vec3D class, 
 * you can find processing here: www.processing.org
 * and you can find the Vec3D class here: code.google.com/p/toxiclibs
 * 
 * @author ruimadeira
 * 
 */
public class MarchingCubes {
	
	PApplet p5;
	
	public float voxelValues[][][];
	protected Vec3D voxels[][][];
	protected Vec3D numPoints, aabbMin, aabbMax;
	protected Vec3D cubeSize;
	protected Vec3D worldSize;
	protected float isoLevel;
	
	private Vec3D vertList[];
	
	protected ArrayList<MCTriangle> triangles;
	
	/**
	 * constructor:
	 * you must define the world bounds, the number of points that will make the grid (in a Vec3D),
	 * and the isoLevel.
	 * @param _p5
	 * @param _aabbMin
	 * @param _aabbMax
	 * @param _numPoints
	 * @param _isoLevel
	 */

	public MarchingCubes(PApplet _p5, Vec3D _aabbMin, Vec3D _aabbMax, Vec3D _numPoints, float _isoLevel){
		p5 = _p5;
		aabbMin = new Vec3D(_aabbMin);
		aabbMax = new Vec3D(_aabbMax);
		worldSize = aabbMax.sub(aabbMin);
		numPoints = new Vec3D(_numPoints);
		cubeSize = new Vec3D(worldSize.x / (numPoints.x-1), worldSize.y / (numPoints.y-1), worldSize.z / (numPoints.z-1));
		voxelValues = new float[(int)numPoints.x][(int)numPoints.y][(int)numPoints.z];
		voxels = new Vec3D[(int)numPoints.x][(int)numPoints.y][(int)numPoints.z];

		_internalReset();
		isoLevel = _isoLevel;
		
		vertList = new Vec3D[12];
		triangles = new ArrayList<MCTriangle>();

	}
	
	/**
	 * creates the mesh
	 */
	public void createMesh(){
		triangles = new ArrayList<MCTriangle>();
		for(int i=0; i<numPoints.x-1; i++){
			for(int j=0; j<numPoints.y-1; j++){
				for(int k=0; k<numPoints.z-1; k++){
					polygonise(i, j, k);
				}
			}
		}
	}
	
	/**
	 * returns an ArrayList of MCTriangles with all the triangles that make up the mesh
	 * @return
	 */
	public ArrayList<MCTriangle> getMesh(){
		return triangles;
	}
	
	/**
	 * copies the mesh triangles into an array and returns
	 * @return
	 */
	public MCTriangle[] getMeshToArray(){
		MCTriangle _triArray[] = new MCTriangle[triangles.size()];
		triangles.toArray(_triArray);
		return _triArray;
	}
	
	/**
	 * default rendering, renders the mesh
	 */
	public void renderMesh(){
		MCTriangle tri;
		p5.beginShape(PApplet.TRIANGLES);
		for(int i=0; i<triangles.size(); i++){
			tri = triangles.get(i);
			p5.vertex(tri.a.x, tri.a.y, tri.a.z);
			p5.vertex(tri.b.x, tri.b.y, tri.b.z);
			p5.vertex(tri.c.x, tri.c.y, tri.c.z);
		}
		p5.endShape();
	}
	
	/**
	 * renders the iso grid.
	 * its useful for debuging.
	 */
	public void renderGrid(){
		p5.noFill();
		p5.stroke(127);
		p5.beginShape(PApplet.LINES);
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z-1; k++){
					p5.vertex(voxels[i][j][k].x, voxels[i][j][k].y, voxels[i][j][k].z);
					p5.vertex(voxels[i][j][k+1].x, voxels[i][j][k+1].y, voxels[i][j][k+1].z);
				}
			}
		}
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y-1; j++){
				for(int k=0; k<numPoints.z; k++){
					p5.vertex(voxels[i][j][k].x, voxels[i][j][k].y, voxels[i][j][k].z);
					p5.vertex(voxels[i][j+1][k].x, voxels[i][j+1][k].y, voxels[i][j+1][k].z);
				}
			}
		}
		
		for(int i=0; i<numPoints.x-1; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z; k++){
					p5.vertex(voxels[i][j][k].x, voxels[i][j][k].y, voxels[i][j][k].z);
					p5.vertex(voxels[i+1][j][k].x, voxels[i+1][j][k].y, voxels[i+1][j][k].z);
				}
			}
		}
		p5.endShape();
	}
	

	/**
	 * returns a tridimensional array of the values that each voxel has.
	 * you can use this to define the value of each voxel
	 * @return
	 */
	public float[][][] getValues(){
		return voxelValues;
	}
	
	/**
	 * return the voxel grid that makes up the iso space, in a three dimensional array.
	 * @return
	 */
	public Vec3D[][][] getVoxels(){
		return voxels;
	}
	
	/**
	 * sets the iso value of a voxel
	 * 
	 * @param posX
	 * @param posY
	 * @param posZ
	 * @param value
	 */
	public void setValue(int indexX, int indexY, int indexZ, float value){
		if(indexX > -1 && indexX < numPoints.x &&
		   indexY > -1 && indexY < numPoints.y && 
		   indexZ > -1 && indexZ < numPoints.z){
			voxelValues[indexX][indexY][indexZ] = value;
		}
	}
	
	/**
	 * gets the value of the specified voxel
	 * @param posX
	 * @param posY
	 * @param posZ
	 * @return
	 */
	public float getValue(int posX, int posY, int posZ){
		if(posX > -1 && posX < numPoints.x &&
		   posY > -1 && posY < numPoints.y && 
		   posZ > -1 && posZ < numPoints.z){
			return voxelValues[posX][posY][posZ];
		} 
		return 0;
	}
	
	/**
	 * returns the a specific voxel of the iso space
	 * @param posX
	 * @param posY
	 * @param posZ
	 * @return
	 */
	public Vec3D getVoxel(int posX, int posY, int posZ){
		if(posX > -1 && posX < numPoints.x &&
		   posY > -1 && posY < numPoints.y && 
		   posZ > -1 && posZ < numPoints.z){
			return voxels[posX][posY][posZ];
		}
		return new Vec3D(0,0,0);
	}
	
	/**
	 * checks if the specified point is inside a voxel cube and returns the voxel.
	 * returns a new Vec3D if point is outside the grid.
	 * @param pos
	 * @return
	 */
	public Vec3D getVoxelAtWorldCoord(Vec3D point){
		for(int i=0; i<voxels.length-1; i++){
			for(int j=0; j<voxels[i].length-1; j++){
				for(int k=0; k<voxels[i][j].length-1; k++){
					if(point.x >= voxels[i][j][k].x && 
					   point.y >= voxels[i][j][k].y && 
					   point.z >= voxels[i][j][k].z &&
					   point.x <= voxels[i+1][j+1][k+1].x &&
					   point.y <= voxels[i+1][j+1][k+1].y &&
					   point.z <= voxels[i+1][j+1][k+1].z){
						return voxels[i][j][k];
					}
				}
			}
		}
		return new Vec3D();
	}
	
	/**
	 *  adds a metaball, with the specified radius, the grid points
	 *  inside the radius will be added the "metaValue"
	 * @param pos
	 * @param radius
	 * @param metaValue
	 */
	public void addMetaBall(Vec3D pos, float radius, float metaValue){
		float radiusSQ = radius*radius;
		float distSQ;
		
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=0; k<voxels[i][j].length; k++){
					distSQ = voxels[i][j][k].distanceToSquared(pos);
					if(distSQ < radiusSQ){
						voxelValues[i][j][k] += (1-distSQ / radiusSQ) * metaValue;
					}
				}
			}
		}
	}
	
	public void addMetaBox(Vec3D aabbMin, Vec3D aabbMax, float metaValue){
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=0; k<voxels[i][j].length; k++){
					if(voxels[i][j][k].x > aabbMin.x && voxels[i][j][k].y > aabbMin.y &&
					   voxels[i][j][k].z > aabbMin.z && voxels[i][j][k].x < aabbMax.x &&
					   voxels[i][j][k].y < aabbMax.y && voxels[i][j][k].z < aabbMax.z){
						PApplet.println("added");
						voxelValues[i][j][k] += metaValue;
					}
				}
			}
		}
	}
	
	/**
	 * returns the maximum voxel value
	 * @return
	 */
	public float getMax(){
		float _max = voxelValues[0][0][0];
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=1; k<voxels[i][j].length; k++){
					if(_max < voxelValues[i][j][k])_max = voxelValues[i][j][k];
				}
			}
		}
		return _max;
	}
	
	/**
	 * returns the lowest voxel value
	 * @return
	 */
	public float getMin(){
		float _min = voxelValues[0][0][0];
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=1; k<voxels[i][j].length; k++){
					if(_min > voxelValues[i][j][k])_min = voxelValues[i][j][k];
				}
			}
		}
		
		return _min;
	}
	
	/**
	 * multiplies all grid values with _val
	 * @param _val
	 */
	public void scale(float _val){
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=0; k<voxels[i][j].length; k++){
					voxelValues[i][j][k] *= _val;
				}
			}
		}
	}
	
	/**
	 * sets all grid values with _val
	 * @param _val
	 */
	public void set(float _val){
		for(int i=0; i<voxels.length; i++){
			for(int j=0; j<voxels[i].length; j++){
				for(int k=0; k<voxels[i][j].length; k++){
					voxelValues[i][j][k] = _val;
				}
			}
		}
	}
	
	/**
	 * sets the grid point with specified index with the value 
	 * @param indexX
	 * @param indexY
	 * @param indexZ
	 * @param val
	 */
	public void set(int indexX, int indexY, int indexZ, float val){
		if(indexX >-1 && indexX < numPoints.x &&
		   indexY >-1 && indexY < numPoints.y &&
		   indexZ >-1 && indexZ < numPoints.z){
			voxelValues[indexX][indexY][indexZ] = val;
		}
	}
	
	/**
	 * normalizes the voxel values
	 */
	public void normalize(){
		float maxVal = 0;
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z; k++){
					if(voxelValues[i][j][k] > maxVal) maxVal = voxelValues[i][j][k];
				}
			}
		}
		float invertMaxVal = 1.0f/maxVal;
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z; k++){
					voxelValues[i][j][k] *= invertMaxVal;
				}
			}
		}
	}
	
	/**
	 * resets the voxel values to zero
	 */
	public void reset(){
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z; k++){
					voxelValues[i][j][k] = 0;
				}
			}
		}
	}
	
	/**
	 * redefines the minimum bounds of the iso space
	 * @param _aabbMin
	 */
	public void setAABBMin(Vec3D _aabbMin){
		aabbMin.set(_aabbMin);
		_internalReset();
	}
	
	/**
	 * returns the minimum bound of the iso space
	 * @return
	 */
	public Vec3D getAABBMin(){
		return aabbMin;
	}
	
	/**
	 * redefines the maximum bound of the iso space
	 * @param _aabbMax
	 */
	public void setAABBMax(Vec3D _aabbMax){
		aabbMax.set(_aabbMax);
		_internalReset();
	}
	
	/**
	 * returns the maximum bound of the iso space
	 * @return
	 */
	public Vec3D getAABBMax(){
		return aabbMax;
	}
	
	/**
	 * returns the number of triangles that make up the mesh
	 * @return
	 */
	public int getNumTriangles(){
		return triangles.size();
	}
	
	/**
	 * returns the iso level
	 * @return
	 */
	public float getIsoLevel(){
		return isoLevel;
	}
	
	/**
	 * sets the iso level
	 * @param _isoLevel
	 */
	public void setIsoLevel(float _isoLevel){
		isoLevel = _isoLevel;
	}
	
	/**
	 * returns the number of vertexes that make up the iso space
	 * in a Vec3D: the x value represents the number of elements along the X axis,
	 * the y value the number of elements along the Y axis and the z value the number
	 * of elements along the Z axis
	 * @return
	 */
	public Vec3D getNumVoxels(){
		return numPoints;
	}
	
	/**
	 * redefines the number of voxels that make up the grid
	 * @param _numPoints
	 */
	public void setNumVoxels(Vec3D _numPoints){
		numPoints.set(_numPoints.x, _numPoints.y, _numPoints.z);
		voxels = new Vec3D[(int)numPoints.x][(int)numPoints.y][(int)numPoints.z];
		voxelValues = new float[(int)numPoints.x][(int)numPoints.y][(int)numPoints.z];
		_internalReset();
	}
	
	/**
	 * returns the size of a single cube of the iso space
	 * @return
	 */
	public Vec3D getCubeSize(){
		return cubeSize;
	}
	
	/**
	 * returns the total size of the iso space
	 * @return
	 */
	public Vec3D getWorldSize(){
		return worldSize;
	}
	
	//Internals
	protected void _internalReset(){
		for(int i=0; i<numPoints.x; i++){
			for(int j=0; j<numPoints.y; j++){
				for(int k=0; k<numPoints.z; k++){
					voxels[i][j][k] = new Vec3D(cubeSize.x * i, cubeSize.y * j, cubeSize.z * k);
					voxels[i][j][k].x += aabbMin.x;
					voxels[i][j][k].y += aabbMin.y;
					voxels[i][j][k].z += aabbMin.z;
					voxelValues[i][j][k] = 0;
					
				}
			}
		}
	}
	
	protected void polygonise(int i, int j, int k){
		int cubeIndex = 0;
		if (voxelValues[i][j][k] < isoLevel) cubeIndex |= 1;
		if (voxelValues[i+1][j][k] < isoLevel) cubeIndex |= 2;
		if (voxelValues[i+1][j+1][k] < isoLevel) cubeIndex |= 4;
		if (voxelValues[i][j+1][k] < isoLevel) cubeIndex |= 8;
		if (voxelValues[i][j][k+1] < isoLevel) cubeIndex |= 16;
		if (voxelValues[i+1][j][k+1] < isoLevel) cubeIndex |= 32;
		if (voxelValues[i+1][j+1][k+1] < isoLevel) cubeIndex |= 64;
		if (voxelValues[i][j+1][k+1] < isoLevel) cubeIndex |= 128;
		/* Cube is entirely in/out of the surface */
		if (MarchingCubesTables.edgeTable[cubeIndex] == 0){
			return;
		}
		

		/* Find the vertices where the surface intersects the cube */
		
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 1) > 0){
			vertList[0] = vertexInterp(isoLevel, voxels[i][j][k], voxels[i+1][j][k], voxelValues[i][j][k] ,voxelValues[i+1][j][k]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 2) > 0){
			vertList[1] = vertexInterp(isoLevel, voxels[i+1][j][k], voxels[i+1][j+1][k], voxelValues[i+1][j][k], voxelValues[i+1][j+1][k]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 4) > 0){
			vertList[2] = vertexInterp(isoLevel, voxels[i+1][j+1][k], voxels[i][j+1][k], voxelValues[i+1][j+1][k], voxelValues[i][j+1][k]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 8) > 0){
			vertList[3] = vertexInterp(isoLevel, voxels[i][j+1][k], voxels[i][j][k], voxelValues[i][j+1][k], voxelValues[i][j][k]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 16) > 0){
			vertList[4] = vertexInterp(isoLevel, voxels[i][j][k+1], voxels[i+1][j][k+1], voxelValues[i][j][k+1], voxelValues[i+1][j][k+1]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 32) > 0){
			vertList[5] = vertexInterp(isoLevel, voxels[i+1][j][k+1], voxels[i+1][j+1][k+1], voxelValues[i+1][j][k+1], voxelValues[i+1][j+1][k+1]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 64) > 0){
			vertList[6] = vertexInterp(isoLevel, voxels[i+1][j+1][k+1], voxels[i][j+1][k+1], voxelValues[i+1][j+1][k+1], voxelValues[i][j+1][k+1]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 128) > 0){
			vertList[7] = vertexInterp(isoLevel, voxels[i][j+1][k+1], voxels[i][j][k+1], voxelValues[i][j+1][k+1], voxelValues[i][j][k+1]); 
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 256) > 0){
			vertList[8] = vertexInterp(isoLevel, voxels[i][j][k], voxels[i][j][k+1], voxelValues[i][j][k], voxelValues[i][j][k+1]);
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 512) > 0){
			vertList[9] = vertexInterp(isoLevel, voxels[i+1][j][k], voxels[i+1][j][k+1], voxelValues[i+1][j][k], voxelValues[i+1][j][k+1]); 
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 1024) > 0){
			vertList[10] = vertexInterp(isoLevel, voxels[i+1][j+1][k], voxels[i+1][j+1][k+1], voxelValues[i+1][j+1][k], voxelValues[i+1][j+1][k+1]); 
		}
		if ((MarchingCubesTables.edgeTable[cubeIndex] & 2048) > 0){
			vertList[11] = vertexInterp(isoLevel,	voxels[i][j+1][k], voxels[i][j+1][k+1], voxelValues[i][j+1][k], voxelValues[i][j+1][k+1]); 
		}
		
		Vec3D vecA;
		Vec3D vecB;
		Vec3D normalVec = new Vec3D();
		for(i=0; MarchingCubesTables.triTable[cubeIndex][i] != -1; i+=3){
			
			vecA = vertList[MarchingCubesTables.triTable[cubeIndex][i+1]].sub(vertList[MarchingCubesTables.triTable[cubeIndex][i]]);
			vecB = vertList[MarchingCubesTables.triTable[cubeIndex][i+2]].sub(vertList[MarchingCubesTables.triTable[cubeIndex][i+1]]);
			normalVec = vecA.cross(vecB);
	
			Vec3D triA = new Vec3D(vertList[MarchingCubesTables.triTable[cubeIndex][i]].x, vertList[MarchingCubesTables.triTable[cubeIndex][i]].y, vertList[MarchingCubesTables.triTable[cubeIndex][i]].z);
			Vec3D triB = new Vec3D(vertList[MarchingCubesTables.triTable[cubeIndex][i+1]].x, vertList[MarchingCubesTables.triTable[cubeIndex][i+1]].y, vertList[MarchingCubesTables.triTable[cubeIndex][i+1]].z);
			Vec3D triC = new Vec3D(vertList[MarchingCubesTables.triTable[cubeIndex][i+2]].x, vertList[MarchingCubesTables.triTable[cubeIndex][i+2]].y, vertList[MarchingCubesTables.triTable[cubeIndex][i+2]].z);
			triangles.add(new MCTriangle(triA, triB, triC, normalVec));
		}
	}
	
	protected Vec3D vertexInterp(float _isoLevel, Vec3D vertice, Vec3D vertice2, float valP1, float valP2){
		 float mu;
		 Vec3D p = new Vec3D();

		   if (Math.abs(isoLevel-valP1) < 0.00001)
		      return(vertice);
		   if (Math.abs(isoLevel-valP2) < 0.00001)
		      return(vertice2);
		   if (Math.abs(valP1-valP2) < 0.00001)
		      return(vertice);
		   mu = (isoLevel - valP1) / (valP2 - valP1);
		   p.x = vertice.x + mu * (vertice2.x - vertice.x);
		   p.y = vertice.y + mu * (vertice2.y - vertice.y);
		   p.z = vertice.z + mu * (vertice2.z - vertice.z);

		   return p;
	}
	
}

