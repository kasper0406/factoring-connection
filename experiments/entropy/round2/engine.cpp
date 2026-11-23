#include "variants.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
int main(int argc, char **argv) {
  if (argc == 2 && std::string(argv[1]) == "names") {
    for (auto &s : names)
      std::cout << s << "\n";
    return 0;
  }
  if (argc != 5)
    return 2;
  std::string mode = argv[1];
  std::vector<int> ids;
  std::istringstream ids_in(argv[2]);
  std::string item;
  while (std::getline(ids_in, item, ','))
    ids.push_back(std::stoi(item));
  std::ifstream f(argv[3], std::ios::binary);
  std::ofstream out(argv[4], std::ios::binary);
  U ab[2];
  std::vector<std::array<U, 2>> input;
  while (f.read((char *)ab, 16))
    input.push_back({ab[0], ab[1]});
  if (!f.eof())
    return 3;
  for (int id : ids) {
    auto begin = std::chrono::steady_clock::now();
    for (auto row : input) {
      W N = W(row[0]) * row[1];
      Trace tr{row[0], row[1]};
      Ring r{N, big(N), tr};
      Out v = variant(id, r);
      uint16_t words[5]{};
      for (unsigned i = 0; i < 4; ++i)
        words[i] = projection(v[i], r);
      if (r.singular)
        std::fill(words, words + 4, 0);
      words[4] = int(tr.split) | (int(r.singular) << 1);
      out.write((char *)words, 10);
      if (mode == "dump") {
        for (auto &w : v)
          std::cout << w << " ";
        std::cout << tr.split << " " << r.singular << "\n";
      }
    }
    if (mode == "time")
      std::cout << id << " "
                << std::chrono::duration<double, std::micro>(
                       std::chrono::steady_clock::now() - begin)
                           .count() /
                       input.size()
                << "\n";
  }
  return out ? 0 : 4;
}
