// qrcode.C — self-contained ROOT macro: URL -> scannable QR code, rendered as a
// TH2D where each dark module is a 2D Gaussian (centre filled, edges fading),
// drawn with COLZ and saved as a PNG.  y-axis "QR Code", x-axis "Scan Me".
//
//   root -l -q 'qrcode.C("https://root.cern")'
//   root -l -q 'qrcode.C("https://example.com/p?q=1","ocean","H",9,0.40,4,"out.png")'
//
// Args: url, scheme, ecc(L|M|Q|H), sub-pixels/module, Gaussian sigma in module
//       units, quiet-zone width (>=4), output file (extension picks the format),
//       x-axis title, y-axis title, optional plot title (empty -> no title).

#include "TH2D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TError.h"

#include <vector>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cctype>

// ===========================================================================
//  QR encoder (byte mode; compact port of Nayuki's public-domain generator)
// ===========================================================================
namespace qr {

// tables indexed [ecc 0..3 = L,M,Q,H][version 1..40]
static const int kEccCwPerBlock[4][41] = {
  {-1, 7,10,15,20,26,18,20,24,30,18,20,24,26,30,22,24,28,30,28,28,28,28,30,30,26,28,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
  {-1,10,16,26,18,24,16,18,22,22,26,30,22,22,24,24,28,28,26,26,26,26,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28},
  {-1,13,22,18,26,18,24,18,22,20,24,28,26,24,20,30,24,28,28,26,30,28,30,30,30,30,28,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
  {-1,17,28,22,16,22,28,26,26,24,28,24,28,22,24,24,30,28,28,26,28,30,24,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30}
};
static const int kNumEccBlocks[4][41] = {
  {-1,1,1,1,1,1,2,2,2,2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8, 8, 9, 9,10,12,12,12,13,14,15,16,17,18,19,19,20,21,22,24,25},
  {-1,1,1,1,2,2,4,4,4,5, 5, 5, 8, 9, 9,10,10,11,13,14,16,17,17,18,20,21,23,25,26,28,29,31,33,35,37,38,40,43,45,47,49},
  {-1,1,1,2,2,4,4,6,6,8, 8, 8,10,12,16,12,17,16,18,21,20,23,23,25,27,29,34,34,35,38,40,43,45,48,51,53,56,59,62,65,68},
  {-1,1,1,2,4,4,4,5,6,8, 8,11,11,16,16,18,16,19,21,25,25,25,34,30,32,35,37,40,42,45,48,51,54,57,60,63,66,70,74,77,81}
};

inline uint8_t rsMul(uint8_t x, uint8_t y) {
  int z = 0;
  for (int i = 7; i >= 0; --i) { z = (z << 1) ^ ((z >> 7) * 0x11D); z ^= ((y >> i) & 1) * x; }
  return (uint8_t)z;
}
inline std::vector<uint8_t> rsDivisor(int degree) {
  std::vector<uint8_t> r(degree, 0); r[degree - 1] = 1; uint8_t root = 1;
  for (int i = 0; i < degree; ++i) {
    for (size_t j = 0; j < r.size(); ++j) { r[j] = rsMul(r[j], root); if (j + 1 < r.size()) r[j] ^= r[j + 1]; }
    root = rsMul(root, 0x02);
  }
  return r;
}
inline std::vector<uint8_t> rsRemainder(const std::vector<uint8_t>& data, const std::vector<uint8_t>& div) {
  std::vector<uint8_t> r(div.size(), 0);
  for (uint8_t b : data) {
    uint8_t factor = b ^ r[0]; r.erase(r.begin()); r.push_back(0);
    for (size_t j = 0; j < r.size(); ++j) r[j] ^= rsMul(div[j], factor);
  }
  return r;
}
inline int numRawDataModules(int ver) {
  int result = (16 * ver + 128) * ver + 64;
  if (ver >= 2) { int n = ver / 7 + 2; result -= (25 * n - 10) * n - 55; if (ver >= 7) result -= 36; }
  return result;
}
inline int numDataCodewords(int ver, int ecl) {
  return numRawDataModules(ver) / 8 - kEccCwPerBlock[ecl][ver] * kNumEccBlocks[ecl][ver];
}
inline void appendBits(unsigned val, int n, std::vector<bool>& bb) {
  for (int i = n - 1; i >= 0; --i) bb.push_back(((val >> i) & 1) != 0);
}
inline bool getBit(long x, int i) { return ((x >> i) & 1) != 0; }

struct Matrix {
  int size = 0, version = 0, ecl = 0;        // ecl: 0=L 1=M 2=Q 3=H
  std::vector<std::vector<bool>> mod, fn;    // [y][x] dark / function-module
  bool at(int x, int y) const { return mod[y][x]; }
  bool ok() const { return size > 0; }
  void setFn(int x, int y, bool dark) {
    if (x < 0 || y < 0 || x >= size || y >= size) return;
    mod[y][x] = dark; fn[y][x] = true;
  }
  std::vector<int> alignPositions() const {
    std::vector<int> res; if (version == 1) return res;
    int n = version / 7 + 2;
    int step = (version == 32) ? 26 : (version * 4 + n * 2 + 1) / (n * 2 - 2) * 2;
    for (int i = 0, pos = size - 7; i < n - 1; ++i, pos -= step) res.insert(res.begin(), pos);
    res.insert(res.begin(), 6); return res;
  }
  void drawFinder(int x, int y) {
    for (int dy = -4; dy <= 4; ++dy) for (int dx = -4; dx <= 4; ++dx) {
      int d = std::max(std::abs(dx), std::abs(dy)), xx = x + dx, yy = y + dy;
      if (xx >= 0 && xx < size && yy >= 0 && yy < size) setFn(xx, yy, d != 2 && d != 4);
    }
  }
  void drawAlign(int x, int y) {
    for (int dy = -2; dy <= 2; ++dy) for (int dx = -2; dx <= 2; ++dx)
      setFn(x + dx, y + dy, std::max(std::abs(dx), std::abs(dy)) != 1);
  }
  void drawFormatBits(int mask) {
    static const int eccFmt[4] = {1, 0, 3, 2};
    int data = eccFmt[ecl] << 3 | mask, rem = data;
    for (int i = 0; i < 10; ++i) rem = (rem << 1) ^ ((rem >> 9) * 0x537);
    int bits = (data << 10 | rem) ^ 0x5412;
    for (int i = 0; i <= 5; ++i) setFn(8, i, getBit(bits, i));
    setFn(8, 7, getBit(bits, 6)); setFn(8, 8, getBit(bits, 7)); setFn(7, 8, getBit(bits, 8));
    for (int i = 9; i < 15; ++i) setFn(14 - i, 8, getBit(bits, i));
    for (int i = 0; i < 8; ++i) setFn(size - 1 - i, 8, getBit(bits, i));
    for (int i = 8; i < 15; ++i) setFn(8, size - 15 + i, getBit(bits, i));
    setFn(8, size - 8, true);
  }
  void drawVersion() {
    if (version < 7) return;
    long rem = version;
    for (int i = 0; i < 12; ++i) rem = (rem << 1) ^ ((rem >> 11) * 0x1F25);
    long bits = (long)version << 12 | rem;
    for (int i = 0; i < 18; ++i) {
      bool bit = getBit(bits, i); int a = size - 11 + i % 3, b = i / 3;
      setFn(a, b, bit); setFn(b, a, bit);
    }
  }
  void drawFunctionPatterns() {
    for (int i = 0; i < size; ++i) { setFn(6, i, i % 2 == 0); setFn(i, 6, i % 2 == 0); }
    drawFinder(3, 3); drawFinder(size - 4, 3); drawFinder(3, size - 4);
    std::vector<int> ap = alignPositions(); int na = (int)ap.size();
    for (int i = 0; i < na; ++i) for (int j = 0; j < na; ++j)
      if (!((i == 0 && j == 0) || (i == 0 && j == na - 1) || (i == na - 1 && j == 0)))
        drawAlign(ap[i], ap[j]);
    drawFormatBits(0); drawVersion();
  }
  void drawCodewords(const std::vector<uint8_t>& data) {
    size_t i = 0;
    for (int right = size - 1; right >= 1; right -= 2) {
      if (right == 6) right = 5;
      for (int vert = 0; vert < size; ++vert) for (int j = 0; j < 2; ++j) {
        int x = right - j; bool up = ((right + 1) & 2) == 0; int y = up ? size - 1 - vert : vert;
        if (!fn[y][x] && i < data.size() * 8) { mod[y][x] = getBit(data[i >> 3], 7 - (int)(i & 7)); ++i; }
      }
    }
  }
  void applyMask(int mask) {
    for (int y = 0; y < size; ++y) for (int x = 0; x < size; ++x) {
      if (fn[y][x]) continue; bool inv = false;
      switch (mask) {
        case 0: inv = (x + y) % 2 == 0; break;
        case 1: inv = y % 2 == 0; break;
        case 2: inv = x % 3 == 0; break;
        case 3: inv = (x + y) % 3 == 0; break;
        case 4: inv = (x / 3 + y / 2) % 2 == 0; break;
        case 5: inv = x * y % 2 + x * y % 3 == 0; break;
        case 6: inv = (x * y % 2 + x * y % 3) % 2 == 0; break;
        case 7: inv = ((x + y) % 2 + x * y % 3) % 2 == 0; break;
      }
      mod[y][x] = mod[y][x] ^ inv;
    }
  }
  long penalty() const {
    long result = 0; const int n = size;
    for (int y = 0; y < n; ++y) { int run = 1;
      for (int x = 1; x < n; ++x) if (mod[y][x] == mod[y][x-1]) { if (++run == 5) result += 3; else if (run > 5) ++result; } else run = 1; }
    for (int x = 0; x < n; ++x) { int run = 1;
      for (int y = 1; y < n; ++y) if (mod[y][x] == mod[y-1][x]) { if (++run == 5) result += 3; else if (run > 5) ++result; } else run = 1; }
    for (int y = 0; y < n - 1; ++y) for (int x = 0; x < n - 1; ++x) {
      bool c = mod[y][x]; if (c == mod[y][x+1] && c == mod[y+1][x] && c == mod[y+1][x+1]) result += 3; }
    for (int y = 0; y < n; ++y) for (int x = 0; x + 7 <= n; ++x) {
      if (!(mod[y][x] && !mod[y][x+1] && mod[y][x+2] && mod[y][x+3] && mod[y][x+4] && !mod[y][x+5] && mod[y][x+6])) continue;
      bool b = true; for (int k = x-4; k < x; ++k) if (k < 0 || mod[y][k]) { b = false; break; }
      bool a = true; for (int k = x+7; k < x+11; ++k) if (k >= n || mod[y][k]) { a = false; break; }
      if (a || b) result += 40; }
    for (int x = 0; x < n; ++x) for (int y = 0; y + 7 <= n; ++y) {
      if (!(mod[y][x] && !mod[y+1][x] && mod[y+2][x] && mod[y+3][x] && mod[y+4][x] && !mod[y+5][x] && mod[y+6][x])) continue;
      bool b = true; for (int k = y-4; k < y; ++k) if (k < 0 || mod[k][x]) { b = false; break; }
      bool a = true; for (int k = y+7; k < y+11; ++k) if (k >= n || mod[k][x]) { a = false; break; }
      if (a || b) result += 40; }
    int dark = 0; for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) if (mod[y][x]) ++dark;
    int total = n * n;
    int k = (int)((std::labs((long)dark * 20 - (long)total * 10) + total - 1) / total) - 1;
    result += (long)k * 10; return result;
  }
};

inline Matrix encode(const std::string& text, int ecl) {
  Matrix q;
  std::vector<uint8_t> bytes(text.begin(), text.end());
  int len = (int)bytes.size();
  int version = -1;
  for (int v = 1; v <= 40; ++v) {
    int cc = (v <= 9) ? 8 : 16;
    if (4 + cc + 8 * len <= numDataCodewords(v, ecl) * 8) { version = v; break; }
  }
  if (version < 0) return q;                 // too long for any version

  std::vector<bool> bb;                       // mode + count + data + pad
  appendBits(0x4, 4, bb);
  appendBits((unsigned)len, (version <= 9) ? 8 : 16, bb);
  for (uint8_t b : bytes) appendBits(b, 8, bb);
  int capBits = numDataCodewords(version, ecl) * 8;
  appendBits(0, std::min(4, capBits - (int)bb.size()), bb);
  appendBits(0, (8 - (int)bb.size() % 8) % 8, bb);
  for (uint8_t pad = 0xEC; (int)bb.size() < capBits; pad ^= (0xEC ^ 0x11)) appendBits(pad, 8, bb);

  std::vector<uint8_t> dataCw(bb.size() / 8, 0);
  for (size_t i = 0; i < bb.size(); ++i) if (bb[i]) dataCw[i >> 3] |= 1 << (7 - (int)(i & 7));

  // Reed-Solomon: split into blocks, append ECC, interleave.
  int numBlocks = kNumEccBlocks[ecl][version], blockEcc = kEccCwPerBlock[ecl][version];
  int rawCw = numRawDataModules(version) / 8;
  int numShort = numBlocks - rawCw % numBlocks, shortLen = rawCw / numBlocks;
  std::vector<std::vector<uint8_t>> blocks;
  std::vector<uint8_t> div = rsDivisor(blockEcc);
  for (int i = 0, k = 0; i < numBlocks; ++i) {
    int dl = shortLen - blockEcc + (i < numShort ? 0 : 1);
    std::vector<uint8_t> dat(dataCw.begin() + k, dataCw.begin() + k + dl); k += dl;
    std::vector<uint8_t> ecc = rsRemainder(dat, div);
    if (i < numShort) dat.push_back(0);
    dat.insert(dat.end(), ecc.begin(), ecc.end());
    blocks.push_back(dat);
  }
  std::vector<uint8_t> allCw;
  for (size_t i = 0; i < blocks[0].size(); ++i) for (size_t j = 0; j < blocks.size(); ++j)
    if (i != (size_t)(shortLen - blockEcc) || j >= (size_t)numShort) allCw.push_back(blocks[j][i]);

  q.version = version; q.ecl = ecl; q.size = 4 * version + 17;
  q.mod.assign(q.size, std::vector<bool>(q.size, false));
  q.fn.assign(q.size, std::vector<bool>(q.size, false));
  q.drawFunctionPatterns();
  q.drawCodewords(allCw);

  int best = 0; long bestPen = -1;             // pick the lowest-penalty mask
  for (int m = 0; m < 8; ++m) {
    q.applyMask(m); q.drawFormatBits(m); long p = q.penalty();
    if (bestPen < 0 || p < bestPen) { bestPen = p; best = m; }
    q.applyMask(m);
  }
  q.applyMask(best); q.drawFormatBits(best);
  return q;
}

// Sub-pixel Gaussian field; v[j*nb+i] with i = x bin (L->R), j = y bin
// (bottom->top, ROOT convention). QR row 0 maps to the top so it is not mirrored.
struct Field { int nb = 0, modules = 0; std::vector<double> v;
  double get(int i, int j) const { return v[(size_t)j * nb + i]; } };

inline Field render(const Matrix& q, int sub, double sigma, int quiet) {
  Field f; f.modules = q.size + 2 * quiet; f.nb = f.modules * sub;
  f.v.assign((size_t)f.nb * f.nb, 0.0);
  const double inv2s2 = 1.0 / (2.0 * sigma * sigma), reach = 3.0 * sigma;
  for (int my = 0; my < q.size; ++my) for (int mx = 0; mx < q.size; ++mx) {
    if (!q.at(mx, my)) continue;
    double cx = (mx + quiet) + 0.5, cy = f.modules - ((my + quiet) + 0.5);
    int iLo = std::max(0, (int)std::floor((cx - reach) * sub));
    int iHi = std::min(f.nb - 1, (int)std::ceil((cx + reach) * sub));
    int jLo = std::max(0, (int)std::floor((cy - reach) * sub));
    int jHi = std::min(f.nb - 1, (int)std::ceil((cy + reach) * sub));
    for (int j = jLo; j <= jHi; ++j) {
      double yy = (j + 0.5) / sub - cy;
      for (int i = iLo; i <= iHi; ++i) {
        double xx = (i + 0.5) / sub - cx;
        f.v[(size_t)j * f.nb + i] += std::exp(-(xx * xx + yy * yy) * inv2s2);
      }
    }
  }
  return f;
}

} // namespace qr

// ===========================================================================
//  Palettes.  `scheme` is either a built-in gradient (value 0 -> background,
//  peak -> dark module), or any ROOT predefined palette by name ("kBird",
//  "kViridis", ...) or number (51..113); a trailing "_r" reverses it.
//  ROOT list: https://root.cern/doc/master/classTColor.html ("High quality...")
// ===========================================================================
struct QrScheme { std::vector<double> s, r, g, b; };

static bool qrCustomScheme(const std::string& n, QrScheme& c) {
  if      (n=="classic")  { c.s={0,1};      c.r={1,0};         c.g={1,0};            c.b={1,0}; }
  else if (n=="inverted") { c.s={0,1};      c.r={0,1};         c.g={0,1};            c.b={0,1}; }
  else if (n=="root")     { c.s={0,1};      c.r={1,0.06};      c.g={1,0.13};         c.b={1,0.42}; }
  else if (n=="matrix")   { c.s={0,1};      c.r={0,0.10};      c.g={0.04,1};         c.b={0,0.18}; }
  else if (n=="ocean")    { c.s={0,1};      c.r={0.88,0.02};   c.g={0.96,0.18};      c.b={1,0.46}; }
  else if (n=="fire")     { c.s={0,0.5,1};  c.r={0,0.75,1};    c.g={0,0.05,0.85};    c.b={0,0.02,0.20}; }
  else if (n=="sunset")   { c.s={0,0.55,1}; c.r={1,0.95,0.30}; c.g={0.98,0.45,0.05}; c.b={0.90,0.20,0.35}; }
  else if (n=="purple")   { c.s={0,1};      c.r={1,0.38};      c.g={1,0.02};         c.b={1,0.55}; }
  else return false;
  return true;
}

static int qrRootPaletteNum(const std::string& n) {
  static const struct { const char* name; int num; } P[] = {
    {"kdeepsea",51},{"kgreyscale",52},{"kdarkbodyradiator",53},{"kblueyellow",54},
    {"krainbow",55},{"kinverteddarkbodyradiator",56},{"kbird",57},{"kcubehelix",58},
    {"kgreenredviolet",59},{"kblueredyellow",60},{"kocean",61},{"kcolorprintableongrey",62},
    {"kalpine",63},{"kaquamarine",64},{"karmy",65},{"katlantic",66},{"kaurora",67},
    {"kavocado",68},{"kbeach",69},{"kblackbody",70},{"kbluegreenyellow",71},{"kbrowncyan",72},
    {"kcmyk",73},{"kcandy",74},{"kcherry",75},{"kcoffee",76},{"kdarkrainbow",77},
    {"kdarkterrain",78},{"kfall",79},{"kfruitpunch",80},{"kfuchsia",81},{"kgreyyellow",82},
    {"kgreenbrownterrain",83},{"kgreenpink",84},{"kisland",85},{"klake",86},
    {"klighttemperature",87},{"klightterrain",88},{"kmint",89},{"kneon",90},{"kpastel",91},
    {"kpearl",92},{"kpigeon",93},{"kplum",94},{"kredblue",95},{"krose",96},{"krust",97},
    {"ksandyterrain",98},{"ksienna",99},{"ksolar",100},{"ksouthwest",101},{"kstarrynight",102},
    {"ksunset",103},{"ktemperaturemap",104},{"kthermometer",105},{"kvalentine",106},
    {"kvisiblespectrum",107},{"kwatermelon",108},{"kcool",109},{"kcopper",110},
    {"kgistearth",111},{"kviridis",112},{"kcividis",113}
  };
  if (!n.empty() && (std::isdigit((unsigned char)n[0]) || n[0]=='-')) return std::atoi(n.c_str());
  for (auto& p : P) if (n == p.name) return p.num;
  return -1;
}

// Set the current palette from `scheme`; return the background colour index.
static int qrApplyPalette(std::string scheme) {
  std::transform(scheme.begin(), scheme.end(), scheme.begin(), ::tolower);
  bool reverse = scheme.size() > 2 && scheme.substr(scheme.size() - 2) == "_r";
  if (reverse) scheme.resize(scheme.size() - 2);

  const int NC = 255;
  QrScheme c;
  if (qrCustomScheme(scheme, c)) {
    TColor::CreateGradientColorTable((int)c.s.size(), c.s.data(), c.r.data(), c.g.data(), c.b.data(), NC);
  } else if (int pal = qrRootPaletteNum(scheme); pal >= 0) {
    gStyle->SetPalette(pal);
    gStyle->SetNumberContours(NC);
  } else {
    ::Warning("qrcode", "unknown scheme '%s'; using classic", scheme.c_str());
    qrCustomScheme("classic", c);
    TColor::CreateGradientColorTable((int)c.s.size(), c.s.data(), c.r.data(), c.g.data(), c.b.data(), NC);
  }
  if (reverse) TColor::InvertPalette();
  return TColor::GetColorPalette(0);
}

// ===========================================================================
//  Macro entry point
// ===========================================================================
void qrcode(const char* url     = "https://root.cern",
            const char* scheme  = "classic",
            const char* ecc     = "H",
            int         sub     = 9,
            double      sigma   = 0.40,
            int         quiet   = 4,
            const char* outfile = "qr.png",
            const char* xtitle  = "Scan Me",
            const char* ytitle  = "QR Code",
            const char* title   = "") {

  std::string e = ecc;
  int ecl = (e=="L"||e=="l") ? 0 : (e=="M"||e=="m") ? 1 : (e=="Q"||e=="q") ? 2 : 3;
  if (quiet < 4) quiet = 4;
  if (sub < 3)   sub   = 3;
  if (sigma <= 0) sigma = 0.40;

  qr::Matrix q = qr::encode(url, ecl);
  if (!q.ok()) { ::Error("qrcode", "text too long for a version-40 QR code"); return; }

  const char* eccNames[4] = {"L", "M", "Q", "H"};
  Printf("QR  : \"%s\"", url);
  Printf("      version %d  (%dx%d modules)  ecc=%s", q.version, q.size, q.size, eccNames[ecl]);

  qr::Field f = qr::render(q, sub, sigma, quiet);

  // Fill every bin (+ a tiny floor) so COLZ paints the whole frame; background
  // bins map to palette colour 0 = the scheme background, never undrawn/white.
  TH2D* h = new TH2D("qr", "", f.nb, 0, f.modules, f.nb, 0, f.modules);
  for (int j = 0; j < f.nb; ++j)
    for (int i = 0; i < f.nb; ++i)
      h->SetBinContent(i + 1, j + 1, f.get(i, j) + 1e-4);

  int bgci = qrApplyPalette(scheme);

  bool hasTitle = title && title[0];
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(hasTitle ? 1 : 0);
  if (hasTitle) h->SetTitle(title);
  h->SetContour(255);
  h->SetMinimum(0.0);
  h->SetMaximum(1.0);                          // a lone dark module saturates to "dark"
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);

  // Margins chosen so the data frame is square (modules stay square) while
  // leaving room for the axis titles and the COLZ palette.
  TCanvas* c = new TCanvas("cqr", "QR code", 950, 900);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.16);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.10);
  c->SetFrameFillColor(bgci);
  h->Draw("COLZ");

  if (outfile && outfile[0]) {
    c->SaveAs(outfile);
    Printf("      saved -> %s", outfile);
  }
}
