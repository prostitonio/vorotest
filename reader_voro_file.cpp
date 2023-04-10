
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;


int main(){
	ifstream fin("voro.vor"); // (ВВЕЛИ НЕ КОРРЕКТНОЕ ИМЯ ФАЙЛА)
	
	char buff[50];
    	if (!fin.is_open()) // если файл не открыт
        	cout << "Файл не может быть открыт!\n"; // сообщить об этом
    	else
    	{	
		int count_mol;
    		fin >> count_mol; // считали первое слово из файла
    		cout << "count_mol "<< count_mol << endl; // считали первое слово из файла

		float vol;
    		fin >> vol; // считали первое слово из файла
    		cout << "vol " << vol << endl; // считали первое слово из файла

	 	int count_nei;
    		fin >> count_nei; 
    		cout << "count_nei " << count_nei << endl; // считали первое слово из файла

		int one_nei;
		vector<int> all_nei_mol;
		cout << "all nei: ";

		for (int i =0; i < 100; i++){ 
			fin >> one_nei;
			cout << one_nei << " ";
			all_nei_mol.push_back(one_nei);
		}
		cout << endl;
		
    	}
	return 0;

}
