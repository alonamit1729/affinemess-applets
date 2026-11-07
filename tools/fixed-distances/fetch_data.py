"""Download the large reproduction-only files that are no longer in the published tree of this repository
(see DATA.md).  They are still stored, byte for byte, in git history; this script fetches them from the pinned
commit below and checks their SHA-256.

    python3 tools/fixed-distances/fetch_data.py            # fetch everything missing
    python3 tools/fixed-distances/fetch_data.py cells/rigid_6_6.json data/strat_4.json   # only some
    python3 tools/fixed-distances/fetch_data.py --list
"""
import hashlib, os, sys, urllib.request
COMMIT = 'ba9fb77dc1f62e81d4d550a6bc5405a674aa97c3'
BASE = f'https://raw.githubusercontent.com/alonamit1729/affinemess-applets/{COMMIT}/tools/fixed-distances/'
FILES = {
    'data/exact_cells2.json': ('7d70a04236846f13b4efd75fdf17986f4bf28ea05f2a1b57f56123f6dc597251', 22884331),
    'data/strat_4.json': ('0d37e81029af6ad487622778e5277553be113ff0a99498aae8bf1be15fe723b4', 14339356),
    'data/strat_5.json': ('f995d94fa944e90884c0ab08b950f2153800aac867739e7b32458a2c5d508caa', 8924440),
    'data/degen_5.json': ('caa05ef4ff65e3f188781f6b36cade73457f3d806c17f1962df9b3d568fdb8b0', 2619957),
    'data/degen_s_5.json': ('e851454a3ebec9d6228aebd9f8aac850b5db6850df619705a19ca663334e1a3d', 1830518),
    'data/degen_s_6.json': ('a22394a29112567bf5716818e54a092ee88e183573d022423bef0586a10fc5f3', 1726047),
    'data/sat56.json': ('dbd653cb2a39ef7f8a99964a251bba8941dcd25c527494be8a84c331dbb53d96', 1489424),
    'cells/cert_6_6.json': ('45cb1ebc5c176e4aa55916d497a77299022e7c8212757cf4a93777a5df671491', 13614925),
    'cells/rigid_6_6.json': ('1d0b51ab6ffa6d726f3aad43ea8214834f0b06cda5f2f13a1433f60089bda1e2', 11732578),
    'cells/cert_7_6.json': ('9db7c3ab6eb17d66fae9a68c877b845fc522767e78e2053b9da42f6351848c42', 2707090),
    'cells/rigid_7_6.json': ('40a426d98f94107494f6e198cd7e297e85312225af0e5aa27a08846af1fee0ac', 2178900),
}
HERE = os.path.dirname(os.path.abspath(__file__))

def fetch(rel):
    sha, size = FILES[rel]; dest = os.path.join(HERE, rel)
    if os.path.exists(dest) and hashlib.sha256(open(dest, 'rb').read()).hexdigest() == sha:
        print('ok      ', rel); return
    print('fetching', rel, '(%.1f MB)' % (size / 1e6), flush=True)
    data = urllib.request.urlopen(BASE + rel, timeout=300).read()
    if hashlib.sha256(data).hexdigest() != sha: sys.exit('checksum mismatch for ' + rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True); open(dest, 'wb').write(data)

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--list' in sys.argv:
        for k, (sha, size) in FILES.items(): print('%-26s %6.1f MB  %s' % (k, size / 1e6, 'present' if os.path.exists(os.path.join(HERE, k)) else 'missing'))
        sys.exit()
    for rel in (args or FILES):
        if rel not in FILES: sys.exit('unknown file ' + rel)
        fetch(rel)
