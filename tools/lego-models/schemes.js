/* Colour schemes for the Lego models, shared by the contents page, the renders, the 3D view and the instructions.
   Every brick has a role: 0 = side A (the whole surface, for one-sided ones), 1 = side B, 2 = rim, 3 = pillar,
   4 = double curve (where an immersed surface crosses itself), 5 and 6 = second and third rims (links).
   A scheme gives each role a Lego colour [hex, name]; roles a model doesn't use are null.
   Images for a scheme live next to the default ones with "-<id>" appended (steps-<id>/, render-<model>-<id>.webp, ...);
   the first scheme of each model is the default and has no suffix.
   A scheme may instead colour every brick on its own: palette [[hex, name], ...], per = one letter per brick
   ("a" = palette[0], ...) in the order of the model's steps file, counts = [[2x4, 2x2], ...] per palette colour.
   style "true": the step pictures and maps show new bricks in their own colours (with bold outlines) rather than
   in the highlight yellow / blue of the default Möbius pages. */
(function () {
  const RED = ["#c91a09", "bright red"], WHITE = ["#f4f4f4", "white"], DBG = ["#6c6e68", "dark bluish grey"];
  const MODELS = {
    mobius: [
      { id: "classic", name: "Red", colors: [RED, null, RED, DBG] },
      { id: "white-rim", name: "White rim", style: "true", colors: [RED, null, WHITE, DBG] },
      { id: "mondrian", name: "Mondrian", style: "true",
        note: "After Mondrian: white with black lines and blocks of red, yellow and blue. One black line sits on the half twist, where the blocks swap sides.", colors: [null, null, null, null],
        palette: [["#f4f4f4", "white"], ["#05131d", "black"], ["#c91a09", "bright red"], ["#f2cd37", "yellow"], ["#0055bf", "blue"]],
        counts: [[1264, 245], [736, 116], [373, 62], [335, 58], [312, 56]],
        per: "bbaaaaaabbbbbbaaaaaaaaaaaaaaabbbbebbaaaaaaaaaaaaaaaabbbbeeeebbaaaaaaaaaaaaaaaaaababbbeeeebbbaaaaaaaaaaaaaaaaabbaabbeeeeeebabaaaaaaaaaaaaaaaaababaabbeeeeeeeeebbaaaaaaaaaaaaaaabbaaaabbeeeeeeeeeebbaabaaabaaaaaaaaaababaaaabbeeeeeeeeeebbbbbbbbbbbbbaaaaabbaaaaabbeeeeeeeeebeebbbbbabaabaabbbbbaaababaaaaaabbeeeeeeeebbabbbaaaaaaaaaabbbbaabbaaaaaaaaabbbeeeeeeeebebaaabbaaaaaaaaaaaaaabbbbabaaaaaaaabbabbeeeebebebaaaaabbaaaaaaaaaaaaaabbbbaaaaaaaaabbaabbbeebeabaaaaaabaaaaaaaaaaaaaaabbbbaaaaaaabbaaaabaebbbeaaaaaaaabbaaaaaaaaaaaaaaaabbabbaaaaaaabbbaaaaababbbbabaaaaaaaabbbbbbbbbaaaaaaaaabababbaaaaabbbbbaaaaaabababaaaaaababbbbbdbbbbbbaaaaaabbaaabbaabbbebbaaaaababbaaaaaaaabbbbbdddddddddbbbbaaababaaabaaabbbeebbaaaaabaabbbbbaabaababaaabdddddddddddbbbaaabbaaaaabbabbeeebabaaaabadbbdababbbabaaaabddddddddddddbbbbaababaaaabbbbeeeeebbaaaabaddbbdabbaabaaaaaaaabdddddddddddddddbbbbbaaaaaaabeeeeeeebaaaaaabbaddddbbbbaaaaaaaaabbddddddddddddddbbbabaaaabbbbeeeeeaaaabadddbdddbbabaaaaaaaaabddddddddddddddddbbbaaabbcbeeeeeeaaaaabaaddbdddbdbbabaaaaaaaaabbdddddddddddddbbbbabbbccbbeeeeeeaaaaaadbddddbbabaabaaaaaadddddddddddbbbbbbcccccbeeeeeeaabaadadddbdadabaabbaaaaadddddddbdbdbbcccccbeeeeeeaaaaabddbbddadabaaaabaadddddbbdddbbbccccbbeeeeeaaaadbbaaddbdaaaaaaababdddddbbbbbcccccbbeeeeaabaadddaadabdabaaaaaadddddbbbbaaabcccccbeeeeeeaadbaadaabddaadaaaaaaddbbabaaabbcccccbeeeeeeeaabaaaddbbddaabaaaaadbabbaaaabbcccbcbeeeeeaaaadbaaddddaabaaaaaabbaabaaaaabcccccbbeeeeaaaaddbabddabaaaaabbaaaaaaabccccceeeeeaaaadddbabdaabaaaabaabaaaabbccccbeeeeeebaaaddaaaabdddaaaabbaaaaaaabbcccbbbeeeeebaaabddbaaaabdaaaaabbaaaaaabccccceeeeebaaaabbdddaaaaadaaabbaaaaaabbcccbbeeeeeddbaaaaaaddaaabaaabbaaaaaaabcccbeeeeeddaaaaaadddbaaaaaabbaaaaaabccccbbeeeeeeddbaaaabbddaaaaaaabbaaaaabcccccbbeeeeebddddbaaaaaaaaaaababbaaaaaaabccccceeeeebddbaaaaaabddaaaaaabbaaaaabccccbeeeeeabddbaaabddbaaaaaabbaaaaaabbccccccbeeeeaabddbaaaddbaaaaaaabbaaaaabccccbeeeeeaaabddddaabaaaababbaaaaaabccccbeeeeeaabdaaadddbaabaabbbbaaaaaabccccbbeeeeeeaaaddaaabdddbaabbbbbaaaaaabccccccbeeeeeaaaaabdaddbadbbbaabbaaaaababccccbeeeeeeaaabdaabddbddbbabaaabbabbaababccccbeeeeeaaabaaaddddbdbaaaaaabbbeabaabbccccbceeeeeaaaababdbbbdabaaaaabbeeeebebbbcccccbeeeeeaaaaaabbbbbabaaaaaaabebeeeeebbcccceeeeeeaaaaaaabbbbaabaaaaaaebeeeeeebbbbcbbeeeeeeaaaabbbaaaababaaaaaeeeeeebaabbbbbeeeeabababbaaabbaaaaaaaaaeeeeeeeebbaaabbbbeeeebbbbcbcaaaabaabaaaaaaaeeeeeeeebbaaaaaabbbbeeebbccccbabaabbaaaaaaaaaeeeeeeebaaaaaabaabbbbccccccabaaaabaabaaaaaaaeeeeeeeebbaaaaaaaababbcccccbaaaaaabbaaaaaaaaaebeeeeeeebbaaaabbaaaababccccccbaaaaaaabbaaaaaaaabbaaaabeeeeeeebbaaaaabbaaaaaacccccbabaaaaaabbaaaaaaabbaaaaaaaaaaabeeeeebbaaaaabaaaaaaccccccbbaaaaaabbaaaaabaaaaaaaaaaaaaaaaaaaebeeebbaaaaaabbaaaaacccccccbaaaaabaabaabaaaaaaaaaaaaaaaaaaeeebebaaaaaaaabbaaaaaacccccbbaaaaaabbabaaaaaaaaaaaaaaaaaabeebbbaaaaaabbaaaaaacccccbcbaaaaaaabbbaaaaaaaaaaaaaaaaaabbbbaaaaaaabaaaaaaccccccbbaaaaaaaabbaaaaaaaaaaaaaaaaaaaaaabbaaaaaaabbaaaaaaacccccccbbbaaaaaabbbabaaaaaaaaaaaabbbbaaaaaaaabaaaaaacccccccbbaaaaabdbdbbbbabaaaaabbbbbbbbaaaaaabbaaaaaaccccccccbbabaabdddddbbbbbbbbbbbbdddbbaaaaabaaaaaaacccccccccbbabbddddddddbddbddddddddddbbaabbbaaaaaaaaccccccccbbbdddddddddddddddddddddbaabbaaaaaaacccccccccbbddddddddddddddddddddddddbabbaaaaaaaaccccccccbbdbddddddddddddddddddddbbabaaaaaaaaccccccbccbbbdddddddddddddddbdbbaaaaaaaaacccccbcccbbbbdbdddddddddbbbbcbbaaaaaaaaccbbccccccbbbbbbbdbbbbbbbbcbccbbaaaaaaacbccccccccccbbbbbbbbbccccccbbaaaabcccccccccccccbcccccccccccbbaaabcccccccccccccccccccccccccccbaaaccccccccccccccccccccccbbaccccccccccccccccccccccbcccccccccccccccccccccccccccccccccccccccccccccccccccc" },
      { id: "ocean", name: "Ocean", style: "true",
        note: "A breaking wave: deep blue along the middle of the band, lighter towards its single edge, where the foam is. The tan pillars are the piers.", colors: [null, null, null, null],
        palette: [["#0a3463", "dark blue"], ["#0055bf", "blue"], ["#078bc9", "dark azure"], ["#36aebf", "medium azure"], ["#008f9b", "dark turquoise"], ["#adc3c0", "light aqua"], ["#f4f4f4", "white"], ["#e4cd9e", "tan"]],
        counts: [[608, 116], [571, 88], [465, 84], [371, 77], [289, 50], [334, 65], [236, 55], [146, 2]],
        per: "hhfffeffhhhhffffffffffffffgfghhhhfffffeefeeeeeffffffhhhhgggfgffeeedededeeeeffggghghhhgggfffeeedddddccddddeefghhgghhgfgegfddddddcccccccccdddefhghgghhggfgfeedddcdcccccccccccdddhhffgghhggfgfdeeddcccdccccccccccccccdhdhefgghhgggeeddcddcbccbcbcbbbbbccccdhhdeffghhgfgegedcdcccbbbbbbbbbbbbbbbbbbcchchddeffghhggefdcddbbbbbabaaabaababbbbbbbhhccdddefgghhggefdfcdcbcbabaaaaaaaaaaaaaaaabbbbhbhbcddefgghhgfgefcdbcbcbaaaaaaaaaaaaaaaaaaaaaahhbbbcddefgghhgfdegccbbabaaaaaaaaaaaaaaaaaaaaaahbhbccddfgghhegdgcebdbbcaaaaaaaaabbabbbbbaaaaaaaaahhaabbcddefgghhgfedecdbcababaaaaaababbbbbbbbbbbbaaaaaahahaabccddefgghhgdffccbababaaaaaabbbbbcbccbcbbbbaaaaahhaaabbcceefghhfgdcdbcabaaaaaaaabcbccccccccbbbbbbbaaahahaaabbcddeefghhfdgcebccaabaaaaaaabcbccdddddcccccccbbbaaaahhaaaaabbcddefghehcdgbbcabaaaaabbaccdcdddddddddccccbbbbbaaahahaaaabbccdeffghhdcdfbbaaaaaaabacbbcdcddeeededdddcccdcccbbbbahhaaaaaaabcddefffhdeffcebbcaaaaababbcedcefeeeeffeeeeeddccccbbhaahaaaaaabbcedffefcdbbaaabaaabacbcdefeffffgffffeeeeedddccbbhhbaaaaaabcceeffebdefbcdaababababadbecdgffgggfffffffffeeddcchbbhaaaaaabbcddeefdfbcdeabaaaabbcdbecdgggggfgfffeeddcchhbbbaaaaabbcdddefdfbcdbbaaabacacbdcfdgfggggggeeddhdhcbbaaaaabccddeffcddebaabbaabacbdbgcegggggfdhhddccbbaaaabbcdeeffcdfbbccdaabbccbecgdgeggggeehhdcbbaaabccddefefcceaabcdbcbacbecgdeggggfehhfddccbbaaaabccddefefbbcfacdbaaccbdfcgeggghfehddcbbaaabbcccdddffdfbccdaabbbadcbedcgfggghhedccbbaaabbbcceefcdefabdfbaabccbddefdfghffhededccbaaabcbddffcdefabbebabcbedcfdhhggfedccbaaabbccdefcdffbabbdbadcbfedfhgfhfedcbbbbaabcccddfbcefaaceccbaabefedhhfgfeedcbbaaabbbccedfbdefbaabdgddbaedcedhhfeedccbaaaabccdefbcdfgcbaabceedcafdehhfeedccbbbaabbcceefaabcdfgccbadfdbedfhhffddcccbaaabbcddfabcdefdcbaabdgfdcfhhfeedcbbbaabbbbccdffaabceedcbaacegffdehhfdcbcbbbaaabbccdffbaaaabcdefdccggedbfhhffdedccbbaabbcceefbaabdffedcbaacefgcfhhfeedcbbaaabcddffcbaabgfecaabcegdeghhegecdcbbaaaababccefdcbaabfecaabcegfedghhgeccbbaaaabbcdfffdcbaaaafdbcdfebghhfgedcbbaaabbbcdffedbafebaaabcebdfgghhgfddbcbaaaabbcddffgfdcaagfdbaaabccfeghhgfgedbbaaaaaabccdfggfdccbafaabbadeggghhgffdcbbbaaaabcbcfeggfdbaecbaabaaccfddgghhgggddbbbaaaabcddgggfcbfdbaaaabacdbfegghhggeecbbaaaaababcdfggfdcbgaaaaaabbcfdgghhgggfddccbaaaaabbcefggeccfdbaaaaabbdcegfghghffecbbaaaaabbcdgggfdcgecbabaaabbcedffghggeeccbaaaaabccedfggedgcbaaaaabcbdcfegfffdccbaaaaabbddgggdgcacbaaaabcccdeefgfffgffedcbbaaaabbcdefggefdbbdaaaabacbccddfffffffedcdbbaaaaaaacbceggggecccbabaabbbccdddfeffeededcbaaaaaabccdefgggecbcabaaaababbccddeefffeeedccbbaaaaabbdeefggecebbaaaaabbbccdcdfeffffefeedddbbaaaabbcdeeggggedgbcbaaaaaaabbcccddeeffffggffeedddccbbbaaaabbccdfgggfedcbabaaaabbbbcccddeeeffffgggggfefeddcdcbbaaaabbcceefggffcccbbabaaaabbccdcdeefffggggggggggfefedddccccbabaaabbccdffggffdcbbbbaaabbbcbcddeeeffggfggggfeedddcccccbbaaaaaabbbcddfgggdeccbbabaaabbbcbccdedeefffffeeddccccbbcbaaaaabbccedfggfeedcbbbaaaaabbbbccddddeeeeededdcccccbbbbaaaaabbcddefgffddccbbbbaaaababbcccdcdddddddddcccbcbbbbbbaaaaaabccddefggfgeecccbbbbaaaaabbbccccccccccccbbbbbbabaaaabbccdeegggeeddccbbabaaaaabbbbcbcbcccbbbbbbbbbbaaaaabcbcddefggeedddccbbbbaaaaaababbbbbbbbbbbbbabaaaaabbbccddeggfefddccccbbbaaaaaaaaaabaabaaaaaaaaaaaabbbbccccdeegfffdedcccbbaaaaaaaaaaaaaaaaaaaaaaaabbcccddeeffffddddcccbbbaaaaaaaaaaaaaaaaaaaaaabbbccccddefffffeddcdccbbaaaaaaaaaaaaaaaaaabbbbcbccdddeeffeeeddcccbbbbaaaaaaaaaaaababbbccccdedefffffeedddcccbbabaaaaaaaaabbbbbcccdcedeffffeeddddccbbbbbabababbabbbbbbcccdddefffffeeedddccccbbbbbbbbbbcbccdcddeeffffeededdcdcccbcbbbbbcccdcdeefffffffeeeecdcdccccccccccdddddefefffgfeeedddddddddddddeeefffgfffffedeededeeeefffffgggfggfffeefeefefgfffgggggggggffggfgffgggggggggggggg" },
    ],
    trefoil: [
      { id: "classic", name: "Orange, azure, white", colors: [["#fe8a18", "orange"], ["#36aebf", "medium azure"], WHITE, DBG] },
      { id: "red-white-black", name: "Red rim, white and black", colors: [WHITE, ["#05131d", "black"], RED, DBG] },
      { id: "red-white-grey", name: "Red rim, white and grey", colors: [WHITE, ["#a0a5a9", "light bluish grey"], RED, DBG] },
    ],
    antioloid: [   // role 2 = the upright circle, role 5 = the level circle
      { id: "classic", name: "Orange, azure", colors: [["#fe8a18", "orange"], ["#36aebf", "medium azure"], WHITE, DBG, null, ["#0a3463", "dark blue"]] },
      { id: "red-white", name: "Red, white and black", colors: [RED, WHITE, ["#05131d", "black"], ["#a0a5a9", "light bluish grey"], null, DBG] },
    ],
    borromean: [   // roles 2, 5 and 6 = the three rings
      { id: "classic", name: "Black, white, grey rings", colors: [["#fe8a18", "orange"], ["#36aebf", "medium azure"], ["#05131d", "black"], ["#a0a5a9", "light bluish grey"], null, WHITE, DBG] },
      { id: "primary", name: "Red, blue, yellow rings", colors: [WHITE, DBG, RED, ["#a0a5a9", "light bluish grey"], null, ["#0055bf", "blue"], ["#f2cd37", "yellow"]] },
    ],
    boy: [   // one-sided: no side B; no pillars; role 4 is the double curve
      { id: "classic", name: "Azure, orange curve", colors: [["#36aebf", "medium azure"], null, WHITE, null, ["#fe8a18", "orange"]] },
      { id: "white-red", name: "White, red rim, black curve", colors: [WHITE, null, RED, null, ["#05131d", "black"]] },
    ],
  };
  function get(model, id) {
    const list = MODELS[model];
    if (id === undefined) { try { id = new URLSearchParams(location.search).get("scheme"); } catch (_) { id = null; } }
    return list.find(s => s.id === id) || list[0];
  }
  const suffix = s => (s.id === "classic" ? "" : "-" + s.id);
  // the URL of this page with ?scheme= set (or removed for the default), keeping the hash
  function url(s, hash) {
    const p = new URLSearchParams(location.search);
    if (s.id === "classic") p.delete("scheme"); else p.set("scheme", s.id);
    const q = p.toString();
    return location.pathname + (q ? "?" + q : "") + (hash === undefined ? location.hash : hash);
  }
  // the colour [hex, name] of brick n (index in the steps file) with role `role`
  const colorOf = (s, n, role) => (s.per ? s.palette[s.per.charCodeAt(n) - 97] : s.colors[role]);
  // distinct colours of a scheme, in role (or palette) order, with the roles that use each
  function distinct(s) {
    if (s.palette) return s.palette.map((c, i) => ({ hex: c[0], name: c[1], roles: [], counts: s.counts[i] }));
    const out = [];
    s.colors.forEach((c, role) => {
      if (!c) return;
      const hit = out.find(o => o.hex === c[0]);
      if (hit) hit.roles.push(role); else out.push({ hex: c[0], name: c[1], roles: [role] });
    });
    return out;
  }
  // a small row of colour dots for a scheme picker
  const dots = s => distinct(s).map(c => `<i style="background:${c.hex}"></i>`).join("");
  // mix a hex colour towards another by t (0..1); light() says whether it's a pale colour that needs a grey edge
  function mix(h, to, t) {
    const a = [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16)), b = [1, 3, 5].map(i => parseInt(to.slice(i, i + 2), 16));
    return "#" + a.map((v, i) => Math.round(v + (b[i] - v) * t).toString(16).padStart(2, "0")).join("");
  }
  const light = h => [1, 3, 5].every(i => parseInt(h.slice(i, i + 2), 16) > 0xd8);
  const edge = h => (light(h) ? "#8a8a86" : mix(h, "#000000", 0.6));
  // pickers are hidden (hidden attribute) when a model has a single scheme; keep that winning over display rules
  try { const st = document.createElement("style"); st.textContent = "[hidden]{display:none !important}"; document.head.appendChild(st); } catch (_) {}
  /* Inside the affinemess.com viewer (an iframe), a link from one of these pages to another would add an entry to
     the browser's joint session history, and the site's "← Applets" button (history.back()) would then step back
     through the pages instead of closing the viewer. So when framed, links within the applet replace the current
     entry instead. Standalone, links behave normally. */
  let framed = false;
  try { framed = window.self !== window.top; } catch (_) { framed = true; }
  if (framed) document.addEventListener("click", e => {
    if (e.defaultPrevented || e.button !== 0 || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
    const a = e.target.closest && e.target.closest("a[href]");
    if (!a || (a.target && a.target !== "_self") || a.hasAttribute("download")) return;
    const u = new URL(a.href, location.href);
    if (u.origin !== location.origin) return;
    e.preventDefault();
    location.replace(u.href);
  });
  window.LEGO_SCHEMES = { MODELS, get, suffix, url, distinct, colorOf, dots, mix, light, edge };
})();
