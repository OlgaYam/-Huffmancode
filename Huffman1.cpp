// Huffman1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cstdio>

using namespace std;

const int size = 100000; //длина выборки
vector<bool> code;                
map<char,vector<bool>> table;  
string alphabet;
map<char,vector<float>> matrix;

map<char, vector<bool>> vec_init;
map<char, map<char, vector<bool>>> tablesForTrees; 


class Node
{
	public:
	 double a;
	 char c;
	 Node *left, *right;
	 
	 Node(){left=right=NULL;}

	 Node(Node *L, Node *R) 
	 {  left =  L;
	    right = R;
	    a = L->a + R->a;  }
};

struct MyCompare
{
    bool operator()(const Node* l, const Node* r) const { return l->a < r->a; }
};
  
Node* init_tree;
map<char, Node * > Trees;

//генерируем текст с помощью марковского источника
void GenerateText(map<char, double> tree)
{
	//находим первый символ по вектору инициализации
	char sym;// = alphabet[0];
	srand(time(NULL));
	double r = (rand()%100)/(100 * 1.0);
	double sum = 0.0;
	for (auto it = tree.begin(); it != tree.end(); ++it)
	{
		sum += it->second;
		if (r < sum)
			sym = it ->first;
	}
	ofstream fileSource("SourceCode.txt", ios::out | ios::binary);
	fileSource << sym; //записываем первый символ в файл
	for (int i = 1; i < size; ++i)
    {
		// генерируем рандомное число
		float r = (rand()%100)/(100 * 1.0);
		//определяем в каком состоянии находимся сейчас
		for (auto it = matrix.begin(); it != matrix.end(); ++it) 
			if (sym == it->first) {//попали сюда -> значит, нашли нужное состояние 
				float m = 0;
				vector<float> pr = it->second;
				//определяем, куда будем переходить дальше
				for (int j = 0; j < pr.size(); ++j) {
					m += pr[j];
					if (r < m) {
						sym = alphabet[j]; 
						break;
					}
				}
				break;
			}
			fileSource << sym; //запись в файл SourceCode
	}
	fileSource.close();
}

void BuildTable(Node *root)
{	
    if (root->left!=NULL) 
                      { code.push_back(0);
                      BuildTable(root->left);}
     
    if (root->right!=NULL)
                       { code.push_back(1);
                       BuildTable(root->right);}
   
    if (root->left==NULL && root->right==NULL) table[root->c]=code;     
    
    code.pop_back();
}

//строим дерево
Node* BuildTree(map<char,double> m)
{
	// записываем начальные узлы в список list     	     
   list<Node*> t;
   for( map<char,double>::iterator itr=m.begin(); itr!=m.end(); ++itr)
   {  
      Node *p = new Node;
      p->c = itr->first;
      p->a = itr->second;  
      t.push_back(p);      
   }	
	
	//  создаем дерево		
	  while (t.size()!=1)
	  {  
		 t.sort(MyCompare());
    
		 Node *SonL = t.front();
		 t.pop_front();
		 Node *SonR = t.front(); 
		 t.pop_front();
     
		 Node *parent = new Node(SonL,SonR); 
		 t.push_front(parent);

	  }
	
	Node *root = t.front();   //root - указатель на вершину дерева
	return root;
}


void GenerationRandom()
{	
	cout << "\nGenerate random" << endl;
	srand(time(NULL));
	vector<float> v;
	for(int i = 0; i != alphabet.length(); ++i)
	{
		double sum = 0.0;
		for(int j = 0; j != alphabet.length()-1; ++j)
		{
			double r = (rand()%100)/(100 * 1.0)*(1-sum);
			sum += r;
			v.push_back(r);
			cout << r << " ";
		}
		double r = 1- sum;
		v.push_back(r);
		cout << r << endl;
		matrix.insert(pair<char, vector<float> > (alphabet[i], v));
		v.clear();
	}
}

//считываем матрицу переходных вероятностей из текстового файла
void ReadFromFile()
{
	float str = 0;
	ifstream matr("matrix.txt", ios::out | ios::binary);
	int i = 0; int k = 0;
	vector<float> v;
	while (!matr.eof())
	{ 
		matr >> str;
		v.push_back(str);
		i++;
		if (i == alphabet.length()) {
			matrix.insert(pair<char, vector<float> > (alphabet[k], v));
			v.clear();
			i = 0;
			k++;
		}
	}
}

//энтропия практическая
void PracticalEntropy(int countBytes, int len)
{
	double entropy2 = double(countBytes*8)/double(len);
    cout << "Practical entropy = " << entropy2 << endl;
}

//Сжимаем по оценке распределения вероятностей
Node* Coder1(string sourceFile, string binaryFile)
{
	//распределение встречающихся символов
	map<char,double> m; 
	//заполняем распределение нулями, так пока символы не встречались
	for(int i = 0; i != alphabet.length(); ++i)
		m.insert ( pair<char,double>(alphabet[i], 0.0) );

	int lenText = 0;
	ifstream fileSource(sourceFile, ios::out | ios::binary);
	string str;
	while (!fileSource.eof())
	{ 
		fileSource >> str; //считали симовл
		for (int i = 0; i != str.length(); ++i)
		{
			++m.find(str[i])->second;
		}
		lenText += str.length();
	}
	for (auto it = m.begin(); it != m.end(); ++it) 
	{
		it->second = it->second/lenText;
	}
	fileSource.close();

	//создали дерево
	Node * root = BuildTree(m);
	// создаем пары 'символ-код'
	table.clear();
   	BuildTable(root); 

	//сжимаем текст
	ifstream fileIn(sourceFile, ios::out | ios::binary);
	ofstream fileBinary(binaryFile, ios::out | ios::binary);
	int count=0; char buf=0;
	int countBytes = 0; int len = 0;
	while (!fileIn.eof())
    {
		fileIn >> str; //считываем из файла SourceCode
		for (int i = 0; i != str.length(); ++i)
		{
			vector<bool> code = table[str[i]];
			//записываем код побитово
			for(int n = 0; n < code.size(); n++)
			{
				buf = buf | code[n]<<(7-count);   
				count++;   
				//если образовался байт, то записываем его в файл
				if (count==8) { 
					count=0;   
					fileBinary<<buf; 
					buf=0; 
					countBytes++;
				} 
			}
		}
		len +=str.length();
	}
	fileIn.close();
	fileBinary.close();

	//Энтропия теоритическая
	double entropy1 = 0;
	for (auto it = m.begin(); it != m.end(); ++it)
	{
		if (it->second != 0.0) {
			double pk = it->second;
			entropy1 += pk*(log(pk)/log(2));
			cout << pk << " ";
		}
	}
	entropy1 *= -1;
	cout << "\nTeoretical entropy = " << entropy1 << endl;

	//практическая энтропия
	PracticalEntropy(countBytes, len);

	return root;
}

void Decoder1(Node *root)
{
	// считывание из файла output.txt и преобразование обратно
	ifstream F("Binary1.txt", ios::in | ios::binary);
	ofstream Out("Decode1.txt", ios::out | ios::binary);
	
	Node *p = root;
	int count=0; char byte; 
	byte = F.get();
	while(!F.eof())
	{   bool b = byte & (1 << (7-count) ) ; 
		if (b) p=p->right; else p=p->left;
		if (p->left==NULL && p->right==NULL) 
		{
			Out<<p->c; 
			p=root;
		}  
		count++;
		if (count==8) {count=0; byte = F.get();}
	}
	
	F.close();
	Out.close();
}

int findFromAlphabet(char c)
{
	for (int i = 0; i != alphabet.length(); ++i)
		if (alphabet[i] == c)
			return i;
	return -1;
}

map<char, double> createVecInit()
{
	//создаем вектор инициализации
	vector<double> pr;
	map<char, double> tree;
	srand(time(NULL));
	double sum = 0.0;
	for(int i = 0; i != alphabet.length()-1; ++i)
	{
		double r = (rand()%100)/(100 * 1.0)*(1-sum);
		sum += r;
		tree.insert(pair<char,double>(alphabet[i], r) );
	}
	double r = 1- sum;
	tree.insert(pair<char,double>(alphabet[alphabet.length()-1], r) );
	init_tree = BuildTree(tree);
	BuildTable(init_tree);
	vec_init = table;
	table.clear();
	return tree;
}


//создание деревьев
void CreateTrees(map<char,vector<float>> m)
{
	for (auto it = m.begin(); it != m.end(); ++it)
	{
		map<char, double> tree;
		for (int i = 0; i != alphabet.length(); ++i)
		{
			tree.insert(pair<char,double>(alphabet[i], it->second[i]) );
		}
		
		Node* node = BuildTree(tree);
		Trees.insert(pair<char, Node*>(it->first, node));
		code.clear();
		BuildTable(node);
		tablesForTrees.insert(pair<char,map<char, vector<bool> > >(it->first, table));
		table.clear();
		tree.clear();
	}
}

void Coding(string sourceFile, string binaryFile)
{
	//сжимаем текст
	ifstream fileIn(sourceFile, ios::out | ios::binary);
	ofstream fileBinary(binaryFile, ios::out | ios::binary);
	int count=0; char buf=0;
	bool flag = true;
	char c1, c2;
	string str;
	int countBytes = 0; int len = 0;
	while (!fileIn.eof())
    {
		fileIn >> str; //считываем из файла SourceCode
		for (int i = 0; i != str.length(); ++i)
		{
			code.clear();
			if (flag) {code = vec_init[str[i]]; c1=str[0]; flag = false;}
			else {
				c2 = str[i];
				//ищем таблицу по предыдущему симовлу
				for (auto it = tablesForTrees.begin(); it != tablesForTrees.end(); ++it)
					if (c1 == it->first)
					{
						//ищем код для текущего символа, учитывая предыдущий
						for (auto jt = it->second.begin(); jt != it->second.end(); ++jt)
						{
							if (c2 == jt->first)
								code = jt->second;
						}
					}
				c1 = c2;
			}
			//записываем код побитово
			for(int n = 0; n < code.size(); n++)
			{
				buf = buf | code[n]<<(7-count);   
				count++;   
				//если образовался байт, то записываем его в файл
				if (count==8) { 
					count=0;   
					fileBinary<<buf; 
					buf=0; 
					countBytes++;
				} 
			}
		}
		len += str.length();
	}
	fileIn.close();
	fileBinary.close();

	//практическая энтропия
	PracticalEntropy(countBytes, len);
}

//Сжимаем, оценивая матрицу переходных вероятностей по выходу источника
void Coder2(string sourceFile, string binaryFile)
{
	//распределение встречающихся символов
	map<char,vector<float> > m;
	vector<float> v;
	for(int i = 0; i != alphabet.length(); ++i) v.push_back(0.0);
	//заполняем распределение нулями, так пока символы не встречались
	for(int i = 0; i != alphabet.length(); ++i)
		m.insert ( pair<char,vector<float>>(alphabet[i], v) );
	int lenText = 0;
	ifstream fileSource("SourceCode.txt", ios::out | ios::binary);
	string str; int i = 1;
	char c1;
	char c2;
	bool flag = true;
	while (!fileSource.eof())
	{ 
		fileSource >> str; //считали симовл
		if (flag) { c1 = str[0]; flag = false;}
		while (i != str.length())
		{
			c2 = str[i];
			int k = findFromAlphabet(c2);
			if (k != -1) ++m.find(c1)->second[k];
			else cout << "error in Coder2" << endl;
			++i;
			c1 = c2;
		}
		i = 0;
	}
	fileSource.close();

	for (auto it = m.begin(); it != m.end(); ++it)
	{
		vector<float> pr = it->second;
		float sum = 0;
		for (int j = 0; j < pr.size(); ++j) {
			sum += pr[j];
		}
		for (int j = 0; j < pr.size(); ++j) {
			pr[j] = pr[j]/sum;
		}
		it->second = pr;
	}

	tablesForTrees.clear();
	table.clear();
	CreateTrees(m);
	Coding(sourceFile, binaryFile);
}

//сжимаем, используя исходную матрицу
void Coder3(string sourceFile, string binaryFile)
{
	tablesForTrees.clear();
	table.clear();
	Trees.clear();
	CreateTrees(matrix);
	Coding(sourceFile, binaryFile);
}


void Decoder(string binaryFile, string decodeFile)
{
	// считывание из файла output.txt и преобразование обратно
	ifstream F(binaryFile, ios::in | ios::binary);
	ofstream Out(decodeFile, ios::out | ios::binary);
	
	Node *p = init_tree;
	int count=0; char byte; 
	byte = F.get();
	while(!F.eof())
	{   
		bool b = byte & (1 << (7-count) ) ; 
		if (b) p=p->right; else p=p->left;
		if (p->left==NULL && p->right==NULL) 
		{
			char ch = p->c;
			Out<<p->c; 
			for (auto it = Trees.begin(); it != Trees.end(); ++it)
				if (ch == it->first) {
					p=it->second;
					break;
				}
		}  
		count++;
		if (count==8) {count=0; byte = F.get();}
	}
	
	F.close();
	Out.close();
}


int main (int argc, char *argv[])
{
		
	cout << "Enter alphabet (without spaces, one line)" << endl;
	cin >> alphabet;

	//alphabet = "abc"; 
	cout << "\nSelect the method of specifying the matrix:" << endl;
	cout << "1 - generate random" << endl;
	cout << "2 - from file matrix.txt" << endl;
	int k;
	cin >> k;
	switch ( k )
	{
	case 1: 
		GenerationRandom(); break;
	case 2:
		ReadFromFile(); break;
	default: 
		cout << "Wrong operation" << endl; 
		system("pause");
		return 0;
	}
  
	ReadFromFile();
	map<char, double> tree = createVecInit();
	GenerateText(tree);
		
	cout << "\nCoding for the evaluation of the probability distribution" << endl;
	Node * root = Coder1("SourceCode.txt", "Binary1.txt");
	Decoder1(root);
	

	cout << "\nCoding, assessing the transition probability matrix for power output" << endl;
	Coder2("SourceCode.txt", "Binary2.txt");
	Decoder("Binary2.txt", "Decode2.txt");

	cout << "\nCoding, using the original matrix" << endl;
	Coder3("SourceCode.txt", "Binary3.txt");
	Decoder("Binary3.txt", "Decode3.txt");
	
	cout << "Decoding is complete" << endl;
	system("pause");
	return 0;
}
