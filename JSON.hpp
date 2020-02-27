#include <iostream>
#include <fstream>
#include <iostream>

#include "Types.hpp"

using namespace std;

void save_object_as_JSON(Alignment* alignment,char* file){
	ofstream myfile;
	myfile.open (file);
	myfile << "{\n";
	myfile << "  \"sequence_1\": {\n";
	myfile << "    \"sequence\": \"" << alignment->sequence_1->sequence << "\",\n";
	myfile << "    \"length\": \"" << alignment->sequence_1->length << "\"\n";
	myfile << "  },\n";
	myfile << "  \"sequence_2\": {\n";
	myfile << "    \"sequence\": \"" << alignment->sequence_1->sequence << "\",\n";
	myfile << "    \"length\": \"" << alignment->sequence_1->length << "\"\n";
	myfile << "  },\n";
	myfile << "  \"parameters\": {\n";
	myfile << "    \"algorithm\": \"" << alignment->parameters->algorithm << "\",\n";
	myfile << "    \"match\": \"" << alignment->parameters->match << "\",\n";
	myfile << "    \"missmatch\": \"" << alignment->parameters->missmatch << "\",\n";
	myfile << "    \"gap\": \"" << alignment->parameters->gap << "\"\n";
	myfile << "  },\n";
	myfile << "  \"parameters\": {\n";
	myfile << "    \"sequence_1\": {\n";
	myfile << "      \"sequence\": \"" << alignment->result->sequence_1->sequence << "\",\n";
	myfile << "      \"length\": \"" << alignment->result->sequence_1->length << "\"\n";
	myfile << "    },\n";
	myfile << "    \"sequence_2\": {\n";
	myfile << "      \"sequence\": \"" << alignment->result->sequence_1->sequence << "\",\n";
	myfile << "      \"length\": \"" << alignment->result->sequence_1->length << "\"\n";
	myfile << "    }\n";
	myfile << "  }\n";
	myfile << "}\n";
}
