#include <bits/stdc++.h>
#include <stdlib.h>
typedef long long int ll;
using namespace std;
#define endl '\n'

int grid_tracking[10000][10000] = {0};

class Node
{
public:
    int id;
    double x;
    double y;
    Node()
    {
    }
    Node(int node_id, double x_pos, double y_pos)
    {
        id = node_id;
        x = x_pos;
        y = y_pos;
    }
};

class Edge
{
public:
    int nodeId1;
    int nodeId2;
    string weight;
    Edge()
    {
    }
    Edge(int node1, int node2, string w)
    {

        nodeId1 = node1;
        nodeId2 = node2;
        weight = w;
    }
};

/* tokenize the string using delimeter ' ' */
vector<string> tokenize(string s, char del)
{
    vector<string> tokens;
    // stringstream class check1
    stringstream check1(s);
    string intermediate;
    // Tokenizing w.r.t. space ','
    while (getline(check1, intermediate, ' '))
    {
        tokens.push_back(intermediate);
    }
    return tokens;
}
class Block
{
    vector<Node> nodes;
    vector<Edge> edges;
    vector<Node> boundary_nodes;
    vector<Edge> boundary_edges;
    Block *overflow;

public:
    Block()
    {
        overflow = NULL;
    }

    void append_overflow_marker_to_file(pair<int, int> cellId);
    void append_node_to_file(Node node, pair<int, int> cellId);
    void append_edge_to_file(Edge edge, pair<int, int> cellId);
    void append_boundary_node_to_file(Node node, pair<int, int> cellId);
    void append_boundary_edge_to_file(Edge edge, pair<int, int> cellId);
    void insert_node_to_block(Node node, pair<int, int> cellId, int blockSize);
    void insert_boundary_node_to_block(Node node, pair<int, int> cellId, int blockSize);
    void insert_boundary_edge_to_block(Edge edge, pair<int, int> cellId, int blockSize);
    void insert_edge_to_block(Edge e, pair<int, int> cellId, int blockSize);
    bool file_exist(string &filename);
    void print_block(pair<int,int> cellId){

        cout << "cellId" << cellId.first << cellId.second << endl;
        Block *node =  this;
        int overflow_count = 0;
        while (node)
        {
            
            if( overflow_count == 0) 
                cout << "--Primary Block => CellId = (" << cellId.first << "," << cellId.second << ").txt" << endl;
            else if (overflow_count > 0)
                cout << "?? Overflow block => CellId = (" << cellId.first << "," << cellId.second << ")_" << overflow_count <<".txt" << endl;
            for (int i = 0; i < node->nodes.size(); i++)
            {
                cout << node->nodes[i].id << " " << node->nodes[i].x <<" "<< node->nodes[i].y << endl;
            }
            cout << "##" << endl;
            for (int i = 0; i < node->edges.size(); i++)
            {
                cout << node->edges[i].nodeId1 << " " << node->edges[i].nodeId2 << " " << node->edges[i].weight << endl;
            }
            cout << "**" << endl;
            for (int i = 0; i < node->boundary_nodes.size(); i++)
            {
                cout << node->boundary_nodes[i].id << " " << node->boundary_nodes[i].x <<" " << node->boundary_nodes[i].y << endl;
            }
            cout << "%%" << endl;
            for (int i = 0; i < node->boundary_edges.size(); i++)
            {
                cout << node->boundary_edges[i].nodeId1 << " " << node->boundary_edges[i].nodeId2 << " " << node->boundary_edges[i].weight << endl;
            }
            node = node->overflow;
            overflow_count++;
        }
    }   
};

class Cell : public Block
{
    pair<int, int> cell_id;
    pair<double, double> cell_min;
    pair<double, double> cell_max;
    int overflow_count;
    bool edge_enabled;
    bool boundary_node_enabled;
    bool boundary_edge_enabled;

public:
    Block *block;
    Cell()
    {

        overflow_count = 0;
        edge_enabled = false;
        boundary_edge_enabled = false;
        boundary_node_enabled = false;
    }
    Cell(pair<int, int> cellId, pair<double, double> cellMin, pair<double, double> cellMax)
    {
        cell_id.first = cellId.first;
        cell_id.second = cellId.second;
        cell_max.first = cellMax.first;
        cell_max.second = cellMax.second;
        cell_min.first = cellMin.first;
        cell_min.second = cellMin.second;
        block = new Block();
        overflow_count = 0;
        edge_enabled = false;
        boundary_node_enabled = false;
        boundary_edge_enabled = false;
    }
    int get_cell_overflow_count()
    {
        return this->overflow_count;
    }
    void update_overflow_count()
    {
        this->overflow_count = this->overflow_count + 1;
    }
    void print_cell()
    {
        if( this ){
            cout << "Cell Id =>";
            print_pair(this->cell_id);
            cout << "Cell Min Cordinates =>"; 
            print_pair(this->cell_min);
            cout << "Cell Max Cordinates =>";
            print_pair(this->cell_max);
            cout << "Print blocks => " << endl;

            if( this->block ){
                this->block->print_block(this->cell_id);
                cout << endl;
            }
        }   
    }

    
    void print_cell_files(){
        cout << "Cell Id =>";
        print_pair(this->cell_id);
        cout << "Cell Min Cordinates =>"; 
        print_pair(this->cell_min);
        cout << "Cell Max Cordinates =>";
        print_pair(this->cell_max);
        cout << "Print blocks => " << endl;
        pair<int,int>cellId = this->cell_id;
        int overflow_count = this->overflow_count;
        string filename  = "(" + to_string(cellId.first) + "," + to_string(cellId.second) + ").txt";
        ifstream infile(filename);
        if( infile.is_open() ){
            string line;
            while (getline(infile, line))
            {
                cout << line << endl;
            }
        }
        infile.close();
        if( overflow_count > 0 ){
            int i = 0;
            while(i <= overflow_count ){
                filename = "(" + to_string(cellId.first) + "," + to_string(cellId.second) + ")_" +to_string(i)+ ".txt";
                ifstream infile(filename);
                if( infile.is_open()){
                    string line;
                    while( getline(infile,line)){
                        cout << line << endl;
                    }
                }
                i++;
                infile.close();
            }
        }
    }
    void print_pair(const pair<int, int> p)
    {
        cout << "<"<< p.first << "," << p.second << ">" << endl;
    }

    bool is_edge_enabled()
    {
        return this->edge_enabled;
    }

    void update_edge_enabled()
    {
        this->edge_enabled = true;
    }

    bool is_boundary_node_enabled()
    {
        return this->boundary_node_enabled;
    }

    void update_boundary_node_enabled()
    {
        this->boundary_node_enabled = true;
    }

    bool is_boundary_edge_enabled()
    {
        return this->boundary_edge_enabled;
    }

    void update_boundary_edge_enabled()
    {
        this->boundary_edge_enabled = true;
    }

    void insert_node_to_cell(Node node, pair<int, int> cellId, int blockSize)
    {
        // cout << "cellId" << endl;
        // print_pair(cellId);
        block->insert_node_to_block(node, cellId, blockSize);
    }

    void insert_boundary_node_to_cell(Node node, pair<int, int> cellId, int blockSize)
    {
        block->insert_boundary_node_to_block(node, cellId, blockSize);
    }
    void insert_boundary_edge_to_cell(Edge e, pair<int, int> cellId, int blockSize)
    {
        block->insert_boundary_edge_to_block(e, cellId, blockSize);
    }

    void insert_edge_to_cell(Edge e, pair<int, int> cellId, int blockSize)
    {

        block->insert_edge_to_block(e, cellId, blockSize);
    }
};

class Grid : public Cell
{
    vector<vector<Cell>> grid;
    map<int, pair<double, double>> node_map;
    double gridMinX, gridMinY, gridMaxX, gridMaxY, gridWidth, gridHeight;
    int cell_size, gridRow, gridCol;
    // bool edge_enabled;
    // map<int, pair<int, int>> boundary_nodes;
    map<pair<int, int>, set<int>> boundary_nodes;
    // set<pair<Node,pair<int,int>>> boundary_nodes;
    // set<Node,pair<int,int>> boundary_nodes;
    // set<tuple<Edge, pair<int, int>, pair<int, int>>> boundary_edges;
    vector<tuple<Edge, pair<int, int>, pair<int, int>>> boundary_edges;
    // bool boundary_node_enabled;
    // bool boundary_edge_enabled;

public:
    Grid()
    {
        // edge_enabled = false;
        // boundary_node_enabled = false;
        // boundary_edge_enabled = false;
    }
    Grid(vector<pair<double, double>> v, int k)
    {
        gridMinX = v[0].first;
        gridMaxX = v[1].first;
        gridMinY = v[0].second;
        gridMaxY = v[1].second;
        gridWidth = gridMaxX - gridMinX;
        gridHeight = gridMaxY - gridMinY;
        gridCol = (gridWidth + k - 1) / k;
        gridRow = (gridHeight + k - 1) / k;
        cell_size = k;
        // edge_enabled = false;
        // boundary_node_enabled = false;
        // boundary_edge_enabled = false;
    }

    void update_overflow_count(pair<int, int> cellId)
    {
        grid[cellId.first][cellId.second].update_overflow_count();
    }

    int get_cell_overflow_count(pair<int, int> cellId)
    {

        return grid[cellId.first][cellId.second].get_cell_overflow_count();
    }

    bool is_edge_enabled(pair<int, int> cellId)
    {

        return grid[cellId.first][cellId.second].is_edge_enabled();
    }

    void update_edge_enabled(pair<int, int> cellId)
    {

        grid[cellId.first][cellId.second].update_edge_enabled();
    }

    bool is_boundary_edge_enabled(pair<int, int> cellId)
    {

        return grid[cellId.first][cellId.second].is_boundary_edge_enabled();
    }

    void update_boundary_edge_enabled(pair<int, int> cellId)
    {

        grid[cellId.first][cellId.second].update_boundary_edge_enabled();
    }

    bool is_boundary_node_enabled(pair<int, int> cellId)
    {

        return grid[cellId.first][cellId.second].is_boundary_node_enabled();
    }

    void update_boundary_node_enabled(pair<int, int> cellId)
    {

        grid[cellId.first][cellId.second].update_boundary_node_enabled();
    }
    void print_grid_size()
    {
        cout << "Total Cells = " << gridRow*gridCol << endl;
        // cout << "gridRow = " << gridRow << endl;
        // cout << "gridCol = " << gridCol << endl;
    }
    void print_grid()
    {
        for (int i = 0; i < grid.size(); i++)
        {
            for (int j = 0; j < grid[i].size(); j++)
            {
                // if( grid[i][j] != NULL )
                    grid[i][j].print_cell();
            }
            cout << endl;
        }
    }

    void initialize_grid(vector<pair<double, double>> v, int k, vector<tuple<int, double, double>> node_tuple)
    {

        gridMinX = v[0].first;
        gridMaxX = v[1].first;
        gridMinY = v[0].second;
        gridMaxY = v[1].second;
        gridWidth = gridMaxX - gridMinX;
        gridHeight = gridMaxY - gridMinY;
        gridCol = (gridWidth + k - 1) / k;
        gridRow = (gridHeight + k - 1) / k;
        cell_size = k;
        double x = gridMinX;
        double y = gridMinY;
        int i = 0;
        int j = 0;

        for (int i = 0; i < node_tuple.size(); i++)
        {

            node_map.insert({get<0>(node_tuple[i]), make_pair(get<1>(node_tuple[i]), get<2>(node_tuple[i]))});
        }

        while (x <= gridMaxX)
        {
            y = gridMinY;
            vector<Cell> v;
            j = 0;
            while (y <= gridMaxY)
            {
                Cell cell(make_pair(i, j), make_pair(x, y), make_pair(x + cell_size, y + cell_size));
                v.push_back(cell);
                y += cell_size;
                j++;
            }
            if (v.size() > 0)
            {
                grid.push_back(v);
            }
            x += cell_size;
            i++;
        }
    }

    pair<int, int> get_cell_id(double x, double y)
    {
        int i, j;
        i = (int)((x - gridMinX) / cell_size);
        j = (int)((y - gridMinY) / cell_size);
        // if( i <= gridRow && j <= gridCol )
        return make_pair(i, j);
    }

    pair<double, double> get_node_cordinates(int nodeId)
    {
        if( node_map.find(nodeId) == node_map.end() ){
            return make_pair(-1,-1);
        }
        else return node_map[nodeId];
    }

    bool invalid_pair(pair<double,double> cord){
        return cord.first < 0 || cord.second < 0;
    }

    bool invalid_cell_id(pair<int,int> cellId){
        return cellId.first < 0 || cellId.second < 0;
    }

    bool compare_pair(pair<int, int> p1, pair<int, int> p2)
    {
        return (p1.first == p2.first && p1.second == p2.second);
    }

    void insert_nodes(string filename, int blockSize)
    {
        string line;
        ifstream datasetfile(filename);
        if (datasetfile.is_open())
        {
            vector<double> x_cordinate;
            vector<double> y_cordinate;
            while (getline(datasetfile, line))
            {
                vector<string> v = tokenize(line, ' ');
                int nodeId = stoi(v[0]);
                double x = stod(v[1]);
                double y = stod(v[2]);
                pair<int, int> cell_id = get_cell_id(x, y);
                // cout << "printing cell id" << endl;
                // print_pair(cell_id);
                Node node(nodeId, x, y);
                grid[cell_id.first][cell_id.second].insert_node_to_cell(node, cell_id, blockSize);
            }
            datasetfile.close();
        }
        else
            cout << "Unable to open file";
    }
    void insert_edges(map<pair<int, int>, string> e_map, int blockSize)
    {
        map<pair<int, int>, string>::iterator itr;
        for (itr = e_map.begin(); itr != e_map.end(); itr++)
        {
            int n1 = itr->first.first;
            int n2 = itr->first.second;
            string e_w = itr->second;
            pair<double, double> node1_cordinate = get_node_cordinates(n1);
            pair<double, double> node2_cordinate = get_node_cordinates(n2);

            pair<int, int> node1_cellId = get_cell_id(node1_cordinate.first, node1_cordinate.second);
            pair<int, int> node2_cellId = get_cell_id(node2_cordinate.first, node2_cordinate.second);

            if (compare_pair(node1_cellId, node2_cellId))
            {

                // both node is in same cell
                // need to add those edges in txt
                Edge e(n1, n2, e_w);
                grid[node1_cellId.first][node1_cellId.second].insert_edge_to_cell(e, node1_cellId, blockSize);
            }
            else
            {
                // Node node1(n1, node1_cordinate.first, node1_cordinate.second);
                // Node node2(n2, node1_cordinate.first, node1_cordinate.second);
                Edge edge(n1, n2, e_w);

                if (boundary_nodes.find(node1_cellId) != boundary_nodes.end())
                {
                    boundary_nodes[node1_cellId].insert(n2);
                }
                else
                {
                    set<int> s;
                    s.insert(n2);
                    boundary_nodes.insert(make_pair(node1_cellId, s));
                }

                if (boundary_nodes.find(node2_cellId) != boundary_nodes.end())
                {
                    boundary_nodes[node2_cellId].insert(n1);
                }
                else
                {
                    set<int> s;
                    s.insert(n1);
                    boundary_nodes.insert(make_pair(node2_cellId, s));
                }

                // boundary_edges.insert(make_tuple(edge, node1_cellId, node2_cellId));
                boundary_edges.push_back(make_tuple(edge, node1_cellId, node2_cellId));
                // grid[node2_cellId.first][node2_cellId.second].insert_boundary_node_to_cell(node1,node2_cellId,blockSize);
                // grid[node1_cellId.first][node1_cellId.second].insert_boundary_node_to_cell(node2,node1_cellId,blockSize);
                // Edge e(n1,n2,v[2]);
                // grid[node1_cellId.first].insert_boundary_edge_to_cell(e,node1_cellId,blockSize);
                // grid[node2_cellId.second].insert_boundary_edge_to_cell(e,node2_cellId,blockSize);
                // insert_node_to_cell()
                // boundary nodes
                // boundary edges
            }
        }
    }

    void insert_boundary_nodes(int blockSize)
    {
        map<pair<int, int>, set<int>>::iterator itr;
        // set<pair<Node,pair<int,int>>>::iterator itr;
        for (itr = boundary_nodes.begin(); itr != boundary_nodes.end(); itr++)
        {
            set<int>::iterator set_itr;
            for (set_itr = itr->second.begin(); set_itr != itr->second.end(); set_itr++)
            {
                pair<double, double> node_cordinate = get_node_cordinates(*set_itr);
                Node node(*set_itr, node_cordinate.first, node_cordinate.second);
                grid[itr->first.first][itr->first.second].insert_boundary_node_to_cell(node, itr->first, blockSize);
            }
        }
    }
    void insert_boundary_edges(int blockSize)
    {
        // set<tuple<Edge,pair<int,int>, pair<int,int>>>::iterator itr;
        // for ( itr = boundary_edges.begin(); itr != boundary_edges.end(); itr++)
        // {
        //     Edge e = get<0>(*itr);
        //     pair<int, int> cellId1 = get<1>(*itr);
        //     pair<int, int> cellId2 = get<2>(*itr);
        //     grid[cellId1.first][cellId1.second].insert_boundary_edge_to_cell(get<0>(*itr), get<1>(*itr), blockSize);
        //     grid[cellId2.first][cellId2.second].insert_boundary_edge_to_cell(get<0>(*itr), get<2>(*itr), blockSize);
        // }

        for (int i = 0; i < boundary_edges.size(); i++)
        {
            Edge e = get<0>(boundary_edges[i]);
            pair<int, int> cellId1 = get<1>(boundary_edges[i]);
            pair<int, int> cellId2 = get<2>(boundary_edges[i]);
            grid[cellId1.first][cellId1.second].insert_boundary_edge_to_cell(get<0>(boundary_edges[i]), get<1>(boundary_edges[i]), blockSize);
            grid[cellId2.first][cellId2.second].insert_boundary_edge_to_cell(get<0>(boundary_edges[i]), get<2>(boundary_edges[i]), blockSize);
        }
    }

    void visualize_grid(){

        cout << "Visulaizing Grid" << endl;
        for(int i = 0; i < grid.size(); i++){
            for( int j = 0; j < grid[i].size(); j++){

                // if( grid[i][j] !=  )
                grid[i][j].print_cell_files();

            }
        }
        cout << "visiualized" << endl;
    }

    void visualize_cell_by_nodeId(int nodeId){
        pair<double,double> node_cords = get_node_cordinates(nodeId);

        if( invalid_pair(node_cords)){
            cout << "invalid node" << endl;    
        }
        else{
            pair<int,int> cellId = get_cell_id(node_cords.first,node_cords.second);
            grid[cellId.first][cellId.second].print_cell_files();
        }
       
    }
};

map<pair<int, int>, string> read_edge_file(string filename)
{
    string line;
    ifstream datasetfile(filename);

    map<pair<int, int>, string> m;
    vector<pair<set<int, int>, string>> vs;
    if (datasetfile.is_open())
    {
        while (getline(datasetfile, line))
        {
            vector<string> v = tokenize(line, ' ');
            int n1 = stoi(v[0]);
            int n2 = stoi(v[1]);
            double edge = stod(v[2]);

            if (m.find(make_pair(n1, n2)) == m.end())
            {
                m.insert({make_pair(n1, n2), v[2]});
            }
        }
        datasetfile.close();
    }
    // else
    //     cout << "Unable to open file";
    return m;
}
vector<tuple<int, double, double>> read_node_file(string filename)
{

    string line;
    ifstream datasetfile(filename);
    int xCells = 0;
    int yCells = 0;
    vector<tuple<int, double, double>> node_tuple;
    map<int, pair<double, double>> node_map;
    if (datasetfile.is_open())
    {
        vector<double> x_cordinate;
        vector<double> y_cordinate;
        while (getline(datasetfile, line))
        {
            vector<string> v = tokenize(line, ' ');
            int nodeId = stoi(v[0]);
            // double x = atof(v[1].c_str());
            // double y = atof(v[2].c_str());
            double x = stod(v[1]);
            double y = stod(v[2]);
            x_cordinate.push_back(x);
            y_cordinate.push_back(y);
            node_tuple.emplace_back(make_tuple(nodeId, x, y));
            node_map.insert({nodeId, make_pair(x, y)});
        }
    }
    return node_tuple;
}

/* read record from dataset.txt file */
vector<pair<double, double>> compute_grid_cordinates(int blockSize, string filename)
{
    string line;
    ifstream datasetfile(filename);
    int xCells = 0;
    int yCells = 0;
    vector<pair<double, double>> v;
    if (datasetfile.is_open())
    {
        vector<double> x_cordinate;
        vector<double> y_cordinate;
        while (getline(datasetfile, line))
        {
            vector<string> v = tokenize(line, ' ');
            int nodeId = stoi(v[0]);
            // double x = atof(v[1].c_str());
            // double y = atof(v[2].c_str());
            double x = stod(v[1]);
            double y = stod(v[2]);
            x_cordinate.push_back(x);
            y_cordinate.push_back(y);
        }
        double maxX = *max_element(x_cordinate.begin(), x_cordinate.end());
        double minX = *min_element(x_cordinate.begin(), x_cordinate.end());
        double maxY = *max_element(y_cordinate.begin(), y_cordinate.end());
        double minY = *min_element(y_cordinate.begin(), y_cordinate.end());
        cout << "maxX = " << maxX << endl;
        cout << "maxY = " << maxY << endl;
        cout << "minX = " << minX << endl;
        cout << "minY = " << minY << endl;
        v.push_back(make_pair(minX, minY));
        v.push_back(make_pair(maxX, maxY));
        datasetfile.close();
    }
    return v;
}

Grid g;

bool Block::file_exist(string &name)
{
    if (FILE *file = fopen(name.c_str(), "r"))
    {
        fclose(file);
        return true;
    }
    else
    {
        return false;
    }
}
void Block::append_boundary_node_to_file(Node node, const pair<int, int> cell_id)
{

    int overflow_count = g.get_cell_overflow_count(cell_id);
    bool edge_enabled = g.is_edge_enabled(cell_id);
    bool boundary_node_enabled = g.is_boundary_node_enabled(cell_id);
    ifstream infile;
    ofstream outfile;
    string filename;
    if (overflow_count > 0)
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "Contents of Overflow Disk block File. => (" + to_string(cell_id.first) +"," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt" << endl;
            outfile << "?? Main Disk Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;

            if( !edge_enabled ){
                outfile << "##" <<endl;
                g.update_edge_enabled(cell_id);
            }
            if (!boundary_node_enabled)
            {
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if( !edge_enabled ){
                outfile << "##" <<endl;
                g.update_edge_enabled(cell_id);
            }
            if (!boundary_node_enabled)
            {
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
        }
        outfile << to_string(node.id) << " " << to_string(node.x) << " " << to_string(node.y) << endl;
        outfile.close();
    }
    else
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ").txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "--Primary Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
            if( !edge_enabled ){
                outfile << "##" <<endl;
                g.update_edge_enabled(cell_id);
            }
            if (!boundary_node_enabled)
            {
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if( !edge_enabled ){
                outfile << "##" <<endl;
                g.update_edge_enabled(cell_id);
            }
            if (!boundary_node_enabled)
            {
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
        }
        outfile << to_string(node.id) << " " << fixed << setprecision(8) << node.x << " " << fixed << setprecision(8) << node.y << endl;
        outfile.close();
    }
}
void Block::append_node_to_file(Node node, const pair<int, int> cell_id)
{
    int overflow_count = g.get_cell_overflow_count(cell_id);
    ifstream infile;
    ofstream outfile;
    string filename;
    if (overflow_count > 0)
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "Content of Overflow Disk block File. => (" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt" << endl;
            outfile << "?? Main Disk Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
        }
        else
        {
            outfile.open(filename, ios_base::app);
        }
        outfile << to_string(node.id) << " " << fixed << setprecision(8) << node.x << " " << fixed << setprecision(8) << node.y << endl;
        outfile.close();
    }
    else
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ").txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "--Primary Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
        }
        else
        {
            outfile.open(filename, ios_base::app);
        }
        outfile << to_string(node.id) << " " << fixed << setprecision(8) << node.x << " " << fixed << setprecision(8) << node.y << endl;
        outfile.close();
    }
}

void Block::append_boundary_edge_to_file(Edge e, const pair<int, int> cell_id)
{
    int overflow_count = g.get_cell_overflow_count(cell_id);
    bool boundary_edge_enabled = g.is_boundary_edge_enabled(cell_id);
    bool edge_enabled = g.is_edge_enabled(cell_id);
    bool boundary_node_enabled = g.is_edge_enabled(cell_id);
    ifstream infile;
    ofstream outfile;
    string filename;
    if (overflow_count > 0)
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "Contents of Overflow Disk block File. => (" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt" << endl;
            outfile << "?? Main Disk Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
            // outfile << "Overflow Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ")_" << overflow_count << ".txt" << endl;

            if( !edge_enabled ){
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }

            if( !boundary_node_enabled ){
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
            if (!boundary_edge_enabled)
            {
                outfile << "%%" << endl;
                g.update_boundary_edge_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if( !edge_enabled ){
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }
            if( !boundary_node_enabled ){
                outfile << "**" << endl;
                g.update_boundary_node_enabled(cell_id);
            }
            if (!boundary_edge_enabled)
            {
                outfile << "%%" << endl;
                g.update_boundary_edge_enabled(cell_id);
            }
        }
        outfile << to_string(e.nodeId1) << " " << to_string(e.nodeId2) << " " << e.weight << endl;
        outfile.close();
    }
    else
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ").txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "--Primary Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
            if (!boundary_edge_enabled)
            {
                outfile << "%%" << endl;
                g.update_boundary_edge_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if (!boundary_edge_enabled)
            {
                outfile << "%%" << endl;
                g.update_boundary_edge_enabled(cell_id);
            }
        }
        //outfile.open(filename, ios_base::app);
        // outfile << std::setprecision(18);
        outfile << to_string(e.nodeId1) << " " << to_string(e.nodeId2) << " " << e.weight << endl;
        outfile.close();
    }
}

void Block::append_edge_to_file(Edge e, const pair<int, int> cell_id)
{
    int overflow_count = g.get_cell_overflow_count(cell_id);
    bool edge_enabled = g.is_edge_enabled(cell_id);
    ifstream infile;
    ofstream outfile;
    string filename;
    if (overflow_count > 0)
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "Contents of Overflow Disk block File. => (" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt" << endl;
            outfile << "?? Main Disk Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
            // outfile << " ?? Overflow Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ")_" << overflow_count << ".txt" << endl;

            if (!edge_enabled)
            {
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if (!edge_enabled)
            {
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }
        }
        outfile << to_string(e.nodeId1) << " " << to_string(e.nodeId2) << " " << e.weight << endl;
        outfile.close();
    }
    else
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ").txt";
        if (!file_exist(filename))
        {
            outfile.open(filename, ios_base::app);
            outfile << "--Primary Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ").txt" << endl;
            if (!edge_enabled)
            {
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }
        }
        else
        {
            outfile.open(filename, ios_base::app);
            if (!edge_enabled)
            {
                outfile << "##" << endl;
                g.update_edge_enabled(cell_id);
            }
        }
        //outfile.open(filename, ios_base::app);
        // outfile << std::setprecision(18);
        outfile << to_string(e.nodeId1) << " " << to_string(e.nodeId2) << " " << e.weight << endl;
        outfile.close();
    }
}

void Block::append_overflow_marker_to_file(pair<int, int> cell_id)
{
    int overflow_count = g.get_cell_overflow_count(cell_id);
    ifstream infile;
    ofstream outfile;
    string filename;
    if (overflow_count > 0)
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ")_" + to_string(overflow_count) + ".txt";
        outfile.open(filename, ios_base::app);
        outfile << "?? Overflow Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ")_" + to_string(overflow_count + 1) + ".txt" << endl;
        outfile.close();
    }
    else
    {
        filename = "(" + to_string(cell_id.first) + "," + to_string(cell_id.second) + ").txt";
        outfile.open(filename, ios_base::app);
        outfile << "?? Overflow Block => CellId = (" << to_string(cell_id.first) << "," << to_string(cell_id.second) << ")_1.txt" << endl;
        outfile.close();
    }
}
void Block::insert_node_to_block(Node node, pair<int, int> cellId, int blockSize)
{
    if (nodes.size() < blockSize)
    {
        nodes.push_back(node);
        append_node_to_file(node, cellId);
    }
    else
    {
        // cout << "overflow" << endl;
        if (overflow == NULL)
        {
            overflow = new Block();
            append_overflow_marker_to_file(cellId);
            grid_tracking[cellId.first][cellId.second]++;
            g.update_overflow_count(cellId);
        }
        overflow->insert_node_to_block(node, cellId, blockSize);
    }
}
void Block::insert_boundary_node_to_block(Node node, pair<int, int> cellId, int blockSize)
{
    if (nodes.size() + edges.size() + boundary_nodes.size() < blockSize)
    {
        boundary_nodes.push_back(node);
        append_boundary_node_to_file(node, cellId);
    }
    else
    {
        if (overflow == NULL)
        {
            overflow = new Block();
            append_overflow_marker_to_file(cellId);
            grid_tracking[cellId.first][cellId.second]++;
            g.update_overflow_count(cellId);
            // append_overflow_marker_to_file(cellId);
            // cell->update_overflow_count();
            // grid[cellId.first][cellId.second].
        }
        overflow->insert_boundary_node_to_block(node, cellId, blockSize);
    }
}

void Block::insert_edge_to_block(Edge e, pair<int, int> cellId, int blockSize)
{
    if (nodes.size() + edges.size() < blockSize)
    {
        edges.push_back(e);
        append_edge_to_file(e, cellId);
    }
    else
    {
        if (overflow == NULL)
        {
            overflow = new Block();
            append_overflow_marker_to_file(cellId);
            g.update_overflow_count(cellId);
            grid_tracking[cellId.first][cellId.second]++;
            // append_overflow_marker_to_file(cellId);
        }
        overflow->insert_edge_to_block(e, cellId, blockSize);
    }
}

void Block::insert_boundary_edge_to_block(Edge e, pair<int, int> cellId, int blockSize)
{
    if (nodes.size() + edges.size() + boundary_nodes.size() + boundary_edges.size() < blockSize)
    {
        boundary_edges.push_back(e);
        append_boundary_edge_to_file(e, cellId);
    }
    else
    {
        if (overflow == NULL)
        {
            overflow = new Block();
            append_overflow_marker_to_file(cellId);

            g.update_overflow_count(cellId);
            grid_tracking[cellId.first][cellId.second]++;
            // append_overflow_marker_to_file(cellId);
        }
        overflow->insert_boundary_edge_to_block(e, cellId, blockSize);
    }
}

int main()
{
    int cellSize, blockSize;
    cout << "Please enter the size of the grid" << endl;
    cin >> cellSize;
    cout << "Please enter the disk block size" << endl;
    cin >> blockSize;
    cout << setprecision(15);
    vector<tuple<int, double, double>> v = read_node_file("nodes.txt");
    vector<pair<double, double>> gridCord = compute_grid_cordinates(blockSize, "nodes.txt");
    g.initialize_grid(gridCord, cellSize, v);
    g.print_grid_size();
    g.insert_nodes("nodes.txt", blockSize);
    map<pair<int, int>, string> edge_map = read_edge_file("edges.txt");
    // map<pair<int, int>, string>::iterator itr;
    // ofstream datasetfile("unique_edges.txt");
    // if (datasetfile.is_open())
    // {
    //     for ( itr = edge_map.begin(); itr != edge_map.end(); itr++)
    //     {
    //         datasetfile << (itr->first).first << " " <<  (itr->first).second << " " <<  itr->second << endl;
    //     }
    //     datasetfile.close();
    // }
    // else
    //     cout << "Unable to open file" << endl;

    cout << "Nodes Inserted" << endl;
    g.insert_edges(edge_map, blockSize);
    cout << "Edges Inserted" << endl;
    g.insert_boundary_nodes(blockSize);
    cout << "Boundary Nodes Inserted" << endl;
    g.insert_boundary_edges(blockSize);
    cout << "Boundary Edge Inserted" << endl;
    int options;
    int nodeId;

    while(1){
        cout << "1. Get cellId from nodeId" << endl;
        cout << "2. Visualize the content of disk block and its associated overflow block for a particular Node." << endl;
        cout << "3. Visulaize the whole grid" << endl;
        cout << "4. Exit Menu" << endl;
        cin >> options;
        if(options == 1){
            cout << "Please Enter NodeId" << endl;
            cin >> nodeId;
            pair<double,double> node_cord = g.get_node_cordinates(nodeId);
            if(!g.invalid_pair(node_cord) ){
                pair<int,int> cellId = g.get_cell_id(node_cord.first,node_cord.second);
                cout << "cellId => " << "<"<< cellId.first << "," << cellId.second << ">" <<endl; 
            }
            else{
                cout << "invalid node" << endl;
            }
        }
        else if( options == 2){
            cout << "Please Enter NodeId" << endl;
            cin >> nodeId;
            g.visualize_cell_by_nodeId(nodeId);
        }
        else if(options == 3){
            g.visualize_grid();
        }
        else if(options == 4){
            cout << "Exit" << endl;
            break;
        }
    }
}
