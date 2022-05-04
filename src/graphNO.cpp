#include "graphNO.hpp"
#include <stdlib.h>
#include <iostream>
#include <string.h>

void GraphNO::initEmptyGraph(Vertex nbVertices) {
   nbVertex=nbVertices;
   matrix.clear();
   matrix.resize(nbVertex *nbVertex);
   matrix.shrink_to_fit();
}

bool GraphNO::isEdge(Vertex v1, Vertex v2){
	return matrix[nbVertex * v1+v2];
}

void GraphNO::addEdge(Vertex v1, Vertex v2){
	 matrix[nbVertex * v1+v2]= true;
	 matrix[nbVertex * v2+v1]= true;
}

void GraphNO::removeEdge(Vertex v1, Vertex v2){
	 matrix[nbVertex * v1+v2]= false;
	 matrix[nbVertex * v2+v1]= false;
}

bool GraphNO::hasNoNeighbor(Vertex v1, const std::vector<Vertex> & v2){
	for(Vertex compt =0; compt < v2.size(); compt++){
		if (isEdge(v1,v2[compt])) return false;
	}
	return true;
}

void GraphNO::markingNeighbors(Vertex v1, std::vector<bool> &v) {
	for(Vertex row =0; row < nbVertex; row++){
		if (!v[row] && isEdge(row,v1)) v[row]=true;
	}
	v[v1]=true;
}


void GraphNO::toComplementGraph() {
   for(Vertex row =0; row < nbVertex; row++){
      for(Vertex col =row+1; col < nbVertex; col++){
        if (isEdge(row,col)) removeEdge(row,col);
        else addEdge(row,col);
      }
      std::cout << std::endl;
   }

}

void GraphNO::display() {
   for(Vertex row =0; row < nbVertex; row++){
      std::cout << "neighbors of " << row << " : ";
      for(Vertex col =0; col < nbVertex; col++){
        if (isEdge(row,col)) std::cout <<  col << " ";
      }
      std::cout << std::endl;
   }
}


Vertex GraphNO::getDegree(Vertex vertex){
   Vertex deg=0;
   for(Vertex row =0; row < nbVertex; row++){
      if (isEdge(row,vertex)) deg++;
   }
   return deg;
}

Vertex GraphNO::getDegreeWithUnmarked(Vertex vertex,	const std::vector<bool> &marked) {
	   Vertex deg=0;
	   for(Vertex row =0; row < nbVertex; row++){
	      if ((!marked[row] )&&isEdge(row,vertex)) deg++;
	   }
	   return deg;

}

vector<Vertex> GraphNO::getNeighbors(Vertex vertex) {
	vector<Vertex> v;
	v.reserve(getNbVertices());
	   for(Vertex row =0; row < nbVertex; row++){
	      if (isEdge(row,vertex)) v.push_back(row);
	   }
	v.shrink_to_fit();
	return v;
}

void GraphNO::intersectNeighbors(Vertex vertex, vector<Vertex> &remaining) {
	vector<Vertex> v;
	v.reserve(remaining.size());

	for (std::vector<Vertex>::iterator it = remaining.begin() ; it != remaining.end(); ++it)
		if (isEdge(vertex,*it)) v.push_back(*it);

	v.shrink_to_fit();
	remaining.resize(v.size());
	remaining=v;
	remaining.shrink_to_fit();

}


void GraphNO::get_params( char * Preamble)
                      /* getting Nr_vert and Nr_edge from the preamble string, 
						 containing Dimacs format "p ??? num num" */
{

	char c, *tmp;
	char * pp = Preamble;
	int stop = 0;
	tmp = (char *)calloc(100, sizeof(char));
	int Nr_vert, Nr_edges;
	Nr_vert = 0;
	
	while (!stop && (c = *pp++) != '\0'){
		switch (c)
		  {
			case 'c':
			  while ((c = *pp++) != '\n' && c != '\0');
			  break;
			  
			case 'p':
			  sscanf(pp, "%s %d %d\n", tmp, &Nr_vert, &Nr_edges);
			  stop = 1;
			  break;
			  
			default:
			  break;
		  }
	}
	
	free(tmp);
	
	std::cout << "import of DIMACS graph with " << Nr_vert << " vertices and " << Nr_edges  << " edges"<< std::endl;
	nbVertex = Nr_vert;

	
}


void GraphNO::importGraphDIMACS( char *file)
{

	int MAX_PREAMBLE = 4000;
	char* Preamble = new char [MAX_PREAMBLE];
	int c, oc;
	char * pp = Preamble;
	int i,j;
	FILE *fp;
	
	if ( (fp=fopen(file,"r"))==NULL )
	  { printf("ERROR: Cannot open infile\n"); exit(10); }

	for(oc = '\0' ;(c = fgetc(fp)) != EOF && (oc != '\n' || c != 'e')
		; oc = *pp++ = c);
 
	ungetc(c, fp); 
	*pp = '\0';
	get_params(Preamble);
	
	initEmptyGraph(nbVertex);

	while ((c = fgetc(fp)) != EOF){
		switch (c)
		  {
			case 'e':
			  if (!fscanf(fp, "%d %d", &i, &j))
				{ printf("ERROR: corrupted inputfile\n"); exit(10); }
			  
			  addEdge(j-1, i-1);
			  break;
			  
			case '\n':
			  
			default:
			  break;
		  }
	}
	/* */
	fclose(fp);
	delete[] Preamble;
	
	std::cout << "end of import of DIMACS graph with " << nbVertex << " vertices"<< std::endl;
	



}


