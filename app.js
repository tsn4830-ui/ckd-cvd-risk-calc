/* ============================================================
   腎臟與心血管風險計算器
   模型係數來源與驗證見 index.html「模型來源與驗證」一節。
   所有計算皆於瀏覽器本機完成。
   ============================================================ */

/* ---------- 單位換算 ---------- */
const CHOL = 38.67;                       // mg/dL -> mmol/L
const mgdl2mmol = v => v / CHOL;
const ngsp2ifcc = p => (p - 2.15) * 10.929;   // HbA1c % -> mmol/mol
const ACR_MGG_PER_MGMMOL = 8.84;          // mg/g = mg/mmol * 8.84

/* ---------- SCORE2 / SCORE2-OP 地區校正因子 ---------- */
const SC = { low:{m:[-0.5699,0.7476], f:[-0.7380,0.7019]},
             mod:{m:[-0.1565,0.8009], f:[-0.3143,0.7701]},
             high:{m:[ 0.3207,0.9360], f:[ 0.5710,0.9369]},
             vh:{m:[ 0.5836,0.8294], f:[ 0.9412,0.8329]} };
const SCOP = { low:{m:[-0.34,1.19], f:[-0.52,1.01]},
               mod:{m:[ 0.01,1.25], f:[-0.10,1.10]},
               high:{m:[ 0.08,1.15], f:[ 0.38,1.09]},
               vh:{m:[ 0.05,0.70], f:[ 0.38,0.69]} };
const REGIONS = [['low','低風險地區'],['mod','中風險地區'],['high','高風險地區'],['vh','極高風險地區']];

/* ---------- SCORE2 / SCORE2-OP ---------- */
function score2(p){                       // p: {sex,age,smoke,sbp,tc,hdl,dm,region}  tc/hdl 單位 mmol/L
  const {sex,age,smoke,sbp,tc,hdl,dm,region} = p;
  let u, s;
  if (age < 70){
    const ca=(age-60)/5, cs=(sbp-120)/20, ct=tc-6, ch=(hdl-1.3)/0.5;
    let lp, S0;
    if (sex==='m'){
      lp = 0.3742*ca + 0.6012*smoke + 0.2777*cs + 0.6457*dm + 0.1458*ct - 0.2698*ch
         - 0.0755*ca*smoke - 0.0255*ca*cs - 0.0281*ca*ct + 0.0426*ca*ch - 0.0983*ca*dm;
      S0 = 0.9605;
    } else {
      lp = 0.4648*ca + 0.7744*smoke + 0.3131*cs + 0.8096*dm + 0.1002*ct - 0.2606*ch
         - 0.1088*ca*smoke - 0.0277*ca*cs - 0.0226*ca*ct + 0.0613*ca*ch - 0.1272*ca*dm;
      S0 = 0.9776;
    }
    u = 1 - Math.pow(S0, Math.exp(lp));
    s = SC[region][sex];
  } else {
    const a = age - 73;
    let lp;
    if (sex==='m'){
      lp = 0.0634*a + 0.4245*dm + 0.3524*smoke + 0.0094*(sbp-150) + 0.0850*(tc-6) - 0.3564*(hdl-1.4)
         - 0.0174*a*dm - 0.0247*a*smoke - 0.0005*a*(sbp-150) + 0.0073*a*(tc-6) + 0.0091*a*(hdl-1.4);
      u = 1 - Math.pow(0.7576, Math.exp(lp - 0.0929));
    } else {
      lp = 0.0789*a + 0.6010*dm + 0.4921*smoke + 0.0102*(sbp-150) + 0.0605*(tc-6) - 0.3040*(hdl-1.4)
         - 0.0107*a*dm - 0.0255*a*smoke - 0.0004*a*(sbp-150) - 0.0009*a*(tc-6) + 0.0154*a*(hdl-1.4);
      u = 1 - Math.pow(0.8082, Math.exp(lp - 0.229));
    }
    s = SCOP[region][sex];
  }
  return 1 - Math.exp(-Math.exp(s[0] + s[1]*Math.log(-Math.log(1-u))));
}

/* ---------- CKD Add-on（CKD-PC）----------
   egfr 必填；acr 可為 null（僅用 eGFR 版本）。
   已與 ckdpcrisk.org 官方計算器 280 組數值逐一比對吻合。            */
function ckdAddon(p){
  const {sex,age,smoke,sbp,tc,hdl,dm,region,egfr,acr} = p;
  const x = score2(p);
  if (egfr == null) return {base:x, addon:null};
  const ca=(age-60)/5, cs=(sbp-120)/20, ct=tc-6, ch=(hdl-1.3)/0.5;
  // 由傳統危險因子推估的「預期 eGFR」
  const ex = 87.8980 - 3.7891*ca - (sex==='f'?0.7023:0) - 0.2941*ct + 1.0960*ch - 0.1364*cs
           + 0.1205*dm + 1.3211*smoke
           + 0.0555*ca*ct + 0.1717*ca*ch + 0.0059*ca*cs - 0.8994*ca*dm + 0.2181*ca*smoke;
  const sa = g => Math.min(g,60)/-15;
  const sb = g => Math.min(Math.max(g-60,0),30)/-15;
  const A1 = sa(egfr)-sa(ex), A2 = sb(egfr)-sb(ex), aa = age-73;
  const d = (age<70) ? (0.4713*A1 + 0.0956*A2 - 0.0802*ca*A1 + 0.0088*ca*A2)
                     : (0.3072*A1 + 0.0942*A2 - 0.0127*aa*A1 - 0.0098*aa*A2);
  let r = 1 - Math.pow(1-x, Math.exp(d));
  if (acr == null) return {base:x, addon:r, acrUsed:false};
  // 預期 log8(ACR)
  const exacr = Math.pow(8,
      1 - 0.0225 + 0.0159*ca + (sex==='f'?0.0308:0) + 0.0185*ct - 0.0274*ch + 0.1339*cs
        + 0.2171*dm + 0.0629*smoke
        - 0.0062*ca*ct + 0.0003*ca*ch + 0.0008*ca*cs - 0.0109*ca*dm + 0.0085*ca*smoke
        + 0.4057*(Math.min(egfr-60,0)/-15) + 0.0597*(Math.min(Math.max(egfr-60,0),30)/-15)
        - 0.0916*(Math.max(egfr-90,0)/-15));
  const b = (age<70) ? 0.2432 : 0.2370;
  const log8 = v => Math.log(v)/Math.log(8);
  r = 1 - Math.pow(1-r, Math.exp(b*(log8(acr)-log8(exacr))));
  return {base:x, addon:r, acrUsed:true};
}

/* ---------- SCORE2-Diabetes ---------- */
function score2dm(p){                     // a1c 單位 mmol/mol；tc/hdl mmol/L
  const {sex,age,smoke,sbp,tc,hdl,dxage,a1c,egfr,region} = p;
  const ca=(age-60)/5, cs=(sbp-120)/20, ct=tc-6, ch=(hdl-1.3)/0.5;
  const cd=(dxage-50)/5, cg=(a1c-31)/9.34, ce=(Math.log(egfr)-4.5)/0.15;
  let lp, S0;
  if (sex==='m'){
    lp = 0.5368*ca + 0.4774*smoke + 0.1322*cs + 0.6457 + 0.1102*ct - 0.1087*ch
       - 0.0672*ca*smoke - 0.0268*ca*cs - 0.0983*ca - 0.0181*ca*ct + 0.0095*ca*ch
       - 0.0998*cd + 0.0955*cg - 0.0591*ce + 0.0058*ce*ce - 0.0134*cg*ca + 0.0115*ce*ca;
    S0 = 0.9605;
  } else {
    lp = 0.6624*ca + 0.6139*smoke + 0.1421*cs + 0.8096 + 0.1127*ct - 0.1568*ch
       - 0.1122*ca*smoke - 0.0167*ca*cs - 0.1272*ca - 0.0200*ca*ct + 0.0186*ca*ch
       - 0.1180*cd + 0.1173*cg - 0.0640*ce + 0.0062*ce*ce - 0.0196*cg*ca + 0.0169*ce*ca;
    S0 = 0.9776;
  }
  const u = 1 - Math.pow(S0, Math.exp(lp));
  const s = SC[region][sex];
  return 1 - Math.exp(-Math.exp(s[0] + s[1]*Math.log(-Math.log(1-u))));
}

/* ---------- KFRE ---------- */
function kfre(p){                         // acr 單位 mg/g
  const {sex,age,egfr,acr,northAmerica,labs} = p;
  const male = sex==='m' ? 1 : 0, la = Math.log(acr);
  const eight = !!labs;
  const lp = eight
    ? (-0.1992*(age/10-7.036) + 0.1602*(male-0.5642) - 0.4919*(egfr/5-7.222) + 0.3364*(la-5.137)
       - 0.3441*(labs.alb-3.997) + 0.2604*(labs.phos-3.916) - 0.07354*(labs.bic-25.57) - 0.2228*(labs.ca-9.355))
    : (-0.2201*(age/10-7.036) + 0.2467*(male-0.5642) - 0.5567*(egfr/5-7.222) + 0.4510*(la-5.137));
  const b2 = eight ? (northAmerica?0.9780:0.9827) : (northAmerica?0.9750:0.9832);
  const b5 = eight ? (northAmerica?0.9301:0.9245) : (northAmerica?0.9240:0.9365);
  return { y2: 1-Math.pow(b2,Math.exp(lp)), y5: 1-Math.pow(b5,Math.exp(lp)) };
}

/* ---------- CKD-EPI 2021（去種族） ---------- */
function ckdepi2021(cr, age, sex){
  const k = sex==='f' ? 0.7 : 0.9, a = sex==='f' ? -0.241 : -0.302;
  return 142 * Math.pow(Math.min(cr/k,1),a) * Math.pow(Math.max(cr/k,1),-1.200)
             * Math.pow(0.9938,age) * (sex==='f'?1.012:1);
}

/* ============================================================
   下拉選單資料
   ============================================================ */
const range=(a,b,s)=>{const o=[];for(let v=a;v<=b+1e-9;v+=s)o.push(+v.toFixed(4));return o;};
const ACR_LIST  = [5,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500];
const KACR_LIST = [3,5,8,10,15,20,25,30,40,50,60,80,100,150,200,250,300,400,500,600,800,1000,1250,1500,2000,2500,3000,4000,5000];
const TC_LIST   = range(100,350,5);       // mg/dL
const HDL_LIST  = range(20,100,1);        // mg/dL
const SBP_LIST  = range(90,200,2);
const A1C_LIST  = range(5.0,15.0,0.1);    // %

const acrBand = v => v<30 ? 'A1' : (v<300 ? 'A2' : 'A3');
const gBand = e => e>=90?'G1':e>=60?'G2':e>=45?'G3a':e>=30?'G3b':e>=15?'G4':'G5';

function fill(el, items, sel){
  el.innerHTML = items.map(([v,t]) => `<option value="${v}"${v==sel?' selected':''}>${t}</option>`).join('');
}
function opts(list, fmt){ return list.map(v=>[v, fmt?fmt(v):String(v)]); }

/* 依單位重新標示（value 不變） */
function relabel(el, list, fmt){
  const cur = el.value;
  fill(el, opts(list, fmt), cur);
  el.value = cur;
}

/* ============================================================
   共用結果元件
   ============================================================ */
function pct(x){ return (x*100).toFixed(1); }
function clsSCORE2(r, age){
  const t = age<50 ? [2.5,7.5] : (age<70 ? [5,10] : [7.5,15]);
  const v = r*100;
  if (v < t[0]) return ['low','低至中度風險', t];
  if (v < t[1]) return ['high','高風險', t];
  return ['vhigh','極高風險', t];
}
function clsSCORE2DM(r){
  const v = r*100;
  if (v<5)  return ['low','低風險'];
  if (v<10) return ['mod','中度風險'];
  if (v<20) return ['high','高風險'];
  return ['vhigh','極高風險'];
}
function bigBox(k, v, cls, x, badge){
  return `<div class="big"><div class="k">${k}</div><div class="v ${cls}">${v}</div>
          <div class="x">${x||''}</div>${badge?`<div class="badge ${cls==='low'?'':cls}">${badge}</div>`:''}</div>`;
}

/* ============================================================
   KFRE 面板
   ============================================================ */
function initKFRE(){
  const $ = id => document.getElementById(id);
  fill($('k_age'),  opts(range(18,100,1)), 65);
  fill($('k_egfr'), opts(range(5,59,1), v=>`${v}（${gBand(v)}）`), 30);
  fill($('k_acr'),  opts(KACR_LIST, v=>`${v}（${acrBand(v)}）`), 300);
  fill($('k_alb'),  opts(range(1.5,5.5,0.1), v=>v.toFixed(1)), 4.0);
  fill($('k_phos'), opts(range(1.5,10.0,0.1), v=>v.toFixed(1)), 3.9);
  fill($('k_bic'),  opts(range(10,35,1)), 24);
  fill($('k_ca'),   opts(range(6.0,12.0,0.1), v=>v.toFixed(1)), 9.2);

  const acrFmt = () => $('k_acrUnit').value==='mgg'
      ? (v=>`${v}（${acrBand(v)}）`)
      : (v=>`${(v/ACR_MGG_PER_MGMMOL).toFixed(1)}（${acrBand(v)}）`);
  $('k_acrUnit').addEventListener('change', ()=>{ relabel($('k_acr'), KACR_LIST, acrFmt()); calc(); });

  function calc(){
    const sex = $('k_sex').value, age = +$('k_age').value,
          egfr = +$('k_egfr').value, acr = +$('k_acr').value,
          na = $('k_region').value === 'na';
    const r4 = kfre({sex,age,egfr,acr,northAmerica:na});
    const labs = { alb:+$('k_alb').value, phos:+$('k_phos').value, bic:+$('k_bic').value, ca:+$('k_ca').value };
    const r8 = kfre({sex,age,egfr,acr,northAmerica:na,labs});
    const use8 = $('k_use8').checked;
    const r = use8 ? r8 : r4;
    const row8 = use8
      ? `<tr class="sel"><td>8 變數（＋白蛋白、磷、HCO₃⁻、鈣）</td><td>${pct(r8.y2)}%</td><td>${pct(r8.y5)}%</td></tr>`
      : `<tr><td>8 變數（＋白蛋白、磷、HCO₃⁻、鈣）</td><td colspan="2" style="text-align:left;color:var(--muted)">展開下方「選填」區、填入四項生化值並勾選後顯示</td></tr>`;

    const c2 = r.y2*100>40 ? 'vhigh' : (r.y2*100>10 ? 'high' : 'low');
    const c5 = r.y5*100>=5 ? 'high' : (r.y5*100>=3 ? 'mod' : 'low');

    let act = [];
    if (r.y5*100 >= 3)  act.push('5 年風險已達 <strong>3–5%</strong> 門檻 → 可考慮<strong>轉介／持續腎臟科追蹤</strong>。');
    if (r.y2*100 > 10)  act.push('2 年風險 <strong>&gt; 10%</strong> → 建議進入<strong>多專科團隊照護</strong>。');
    if (r.y2*100 > 40)  act.push('2 年風險 <strong>&gt; 40%</strong> → 應開始<strong>腎臟替代治療衛教、血管通路規劃或轉介移植評估</strong>。');
    if (!act.length)    act.push('目前未達 KDIGO 建議的轉介或準備透析門檻，仍應依 eGFR／白蛋白尿分級持續追蹤與控制危險因子。');

    $('k_result').innerHTML = `
      <div class="big-row">
        ${bigBox('2 年腎衰竭風險', pct(r.y2)+'%', c2, `${use8?'8':'4'} 變數版本・${na?'北美':'非北美'}校正`)}
        ${bigBox('5 年腎衰竭風險', pct(r.y5)+'%', c5, `CKD ${gBand(egfr)}／${acrBand(acr)}`)}
      </div>
      <div class="interpret">${act.map(a=>'• '+a).join('<br>')}</div>
      <table class="regions">
        <thead><tr><th>版本</th><th>2 年</th><th>5 年</th></tr></thead>
        <tbody>
          <tr class="${use8?'':'sel'}"><td>4 變數（年齡、性別、eGFR、uACR）</td><td>${pct(r4.y2)}%</td><td>${pct(r4.y5)}%</td></tr>
          ${row8}
        </tbody>
      </table>
      ${egfr>=60?'<div class="warn">KFRE 是在 CKD G3–G5 族群開發的，eGFR ≥ 60 時結果不可靠。</div>':''}`;
  }
  // 8 變數開關（放進 optional 區塊）
  const opt = document.querySelector('#panel-kfre details.optional');
  const sw = document.createElement('label');
  sw.className = 'optin';
  sw.innerHTML = '<input type="checkbox" id="k_use8"> 我已填好下方四項數值，改用 <strong>8 變數版本</strong>作為主要結果';
  opt.insertBefore(sw, opt.querySelector('.grid'));

  document.querySelectorAll('#panel-kfre select, #panel-kfre input').forEach(el=>el.addEventListener('change', calc));
  calc();
}

/* ============================================================
   SCORE2 / SCORE2-OP / SCORE2-Diabetes 面板
   ============================================================ */
function buildFields(host, pfx, kind){
  const ageRange = kind==='op' ? [70,89] : [40,69];
  const rows = [];
  rows.push(`<div class="field"><label>年齡</label><select id="${pfx}_age"></select><span class="unit">歲</span></div>`);
  rows.push(`<div class="field"><label>性別</label><select id="${pfx}_sex"><option value="m">男</option><option value="f">女</option></select></div>`);
  rows.push(`<div class="field"><label>目前吸菸</label><select id="${pfx}_smoke"><option value="0">否</option><option value="1">是</option></select></div>`);
  rows.push(`<div class="field"><label>收縮壓 SBP</label><select id="${pfx}_sbp"></select><span class="unit">mmHg</span></div>`);
  rows.push(`<div class="field"><label>總膽固醇 TC</label><select id="${pfx}_tc"></select>
             <select id="${pfx}_cholUnit" class="unit-sel"><option value="mgdl">mg/dL</option><option value="mmol">mmol/L</option></select></div>`);
  rows.push(`<div class="field"><label>高密度脂蛋白 HDL-C</label><select id="${pfx}_hdl"></select><span class="unit">單位同上</span></div>`);
  if (kind==='dm'){
    rows.push(`<div class="field"><label>糖尿病診斷年齡</label><select id="${pfx}_dxage"></select><span class="unit">歲</span></div>`);
    rows.push(`<div class="field"><label>HbA1c</label><select id="${pfx}_a1c"></select>
               <select id="${pfx}_a1cUnit" class="unit-sel"><option value="ngsp">%（NGSP）</option><option value="ifcc">mmol/mol（IFCC）</option></select></div>`);
    rows.push(`<div class="field"><label>eGFR</label><select id="${pfx}_egfr"></select><span class="unit">mL/min/1.73m²（必填）</span></div>`);
  } else {
    rows.push(`<div class="field"><label>糖尿病</label><select id="${pfx}_dm"><option value="0">無</option><option value="1">有</option></select></div>`);
    rows.push(`<div class="field"><label>eGFR</label><select id="${pfx}_egfr"></select><span class="unit">mL/min/1.73m²</span></div>`);
    rows.push(`<div class="field"><label>尿液白蛋白／肌酸酐比 uACR</label><select id="${pfx}_acr"></select>
               <select id="${pfx}_acrUnit" class="unit-sel"><option value="mgg">mg/g</option><option value="mgmmol">mg/mmol</option></select></div>`);
  }
  rows.push(`<div class="field"><label>風險地區（歐洲校正）</label><select id="${pfx}_region">
             ${REGIONS.map(([v,t])=>`<option value="${v}"${v==='mod'?' selected':''}>${t}</option>`).join('')}</select></div>`);
  host.innerHTML = rows.join('');

  const $ = id => document.getElementById(id);
  fill($(`${pfx}_age`), opts(range(ageRange[0],ageRange[1],1)), kind==='op'?75:60);
  fill($(`${pfx}_sbp`), opts(SBP_LIST), 140);
  const cholFmt = pfx2 => () => $(`${pfx}_cholUnit`).value==='mgdl' ? (v=>String(v)) : (v=>mgdl2mmol(v).toFixed(2));
  fill($(`${pfx}_tc`),  opts(TC_LIST),  200);
  fill($(`${pfx}_hdl`), opts(HDL_LIST), 50);
  $(`${pfx}_cholUnit`).addEventListener('change', ()=>{
    const f = cholFmt()();
    relabel($(`${pfx}_tc`), TC_LIST, f);
    relabel($(`${pfx}_hdl`), HDL_LIST, f);
  });
  if (kind==='dm'){
    fill($(`${pfx}_dxage`), opts(range(15,69,1)), 50);
    fill($(`${pfx}_a1c`), opts(A1C_LIST, v=>v.toFixed(1)), 7.5);
    fill($(`${pfx}_egfr`), opts(range(15,140,1), v=>`${v}（${gBand(v)}）`), 90);
    $(`${pfx}_a1cUnit`).addEventListener('change', ()=>{
      const ifcc = $(`${pfx}_a1cUnit`).value==='ifcc';
      relabel($(`${pfx}_a1c`), A1C_LIST, ifcc ? (v=>Math.round(ngsp2ifcc(v))) : (v=>v.toFixed(1)));
    });
  } else {
    fill($(`${pfx}_egfr`), [[ '', '未檢驗（不計 Add-on）' ]].concat(opts(range(15,140,1), v=>`${v}（${gBand(v)}）`)), 60);
    fill($(`${pfx}_acr`),  [[ '', '未檢驗（僅用 eGFR）' ]].concat(opts(ACR_LIST, v=>`${v}（${acrBand(v)}）`)), 30);
    $(`${pfx}_acrUnit`).addEventListener('change', ()=>{
      const mgmmol = $(`${pfx}_acrUnit`).value==='mgmmol';
      const cur = $(`${pfx}_acr`).value;
      fill($(`${pfx}_acr`), [[ '', '未檢驗（僅用 eGFR）' ]].concat(
        opts(ACR_LIST, mgmmol ? (v=>`${(v/ACR_MGG_PER_MGMMOL).toFixed(1)}（${acrBand(v)}）`)
                              : (v=>`${v}（${acrBand(v)}）`))), cur);
      $(`${pfx}_acr`).value = cur;
    });
  }
}

function readCommon(pfx){
  const $ = id => document.getElementById(id);
  return {
    age:+$(`${pfx}_age`).value, sex:$(`${pfx}_sex`).value,
    smoke:+$(`${pfx}_smoke`).value, sbp:+$(`${pfx}_sbp`).value,
    tc: mgdl2mmol(+$(`${pfx}_tc`).value), hdl: mgdl2mmol(+$(`${pfx}_hdl`).value),
    region: $(`${pfx}_region`).value
  };
}

function regionTable(rowsFn, selected){
  const rows = REGIONS.map(([k,label])=>{
    const cells = rowsFn(k);
    return `<tr class="${k===selected?'sel':''}"><td>${label}</td>${cells.map(c=>`<td>${c}</td>`).join('')}</tr>`;
  }).join('');
  return rows;
}

function initSCORE2(pfx, panelId, kind){
  buildFields(document.getElementById(`g_${pfx}`), pfx, kind);
  const $ = id => document.getElementById(id);

  function calc(){
    const c = readCommon(pfx);
    const dm = +$(`${pfx}_dm`).value;
    const egfrRaw = $(`${pfx}_egfr`).value, acrRaw = $(`${pfx}_acr`).value;
    const egfr = egfrRaw==='' ? null : +egfrRaw;
    const acr  = acrRaw==='' ? null : +acrRaw;
    const p = {...c, dm, egfr, acr};
    const res = ckdAddon(p);
    const shown = res.addon==null ? res.base : res.addon;
    const [cls, label, thr] = clsSCORE2(shown, c.age);
    const addonTxt = res.addon==null ? '未加 CKD Add-on'
                   : (res.acrUsed ? 'eGFR + uACR Add-on' : '僅 eGFR Add-on');
    const delta = res.addon==null ? '' :
      `<div class="delta">相較原始 SCORE2 ${res.addon>=res.base?'↑':'↓'} ${Math.abs((res.addon-res.base)*100).toFixed(1)} 個百分點</div>`;

    const mainTitle = res.addon==null
      ? `10 年心血管風險（SCORE2${kind==='op'?'-OP':''}）`
      : '10 年心血管風險（含 CKD Add-on）';
    $(`${pfx}_result`).innerHTML = `
      <div class="big-row">
        ${bigBox(mainTitle, pct(shown)+'%', cls, addonTxt, label)}
        ${res.addon==null ? '' : bigBox('原始 SCORE2'+(kind==='op'?'-OP':''), pct(res.base)+'%', 'low', '未考慮腎功能')}
      </div>
      ${delta}
      <div class="interpret">
        <strong>${c.age} 歲的分級門檻</strong>：低至中度 &lt; ${thr[0]}%、高風險 ${thr[0]}–&lt;${thr[1]}%、極高風險 ≥ ${thr[1]}%（ESC 2021）。
        ${dm?'<br>此人有糖尿病 → 建議改用 <strong>SCORE2-Diabetes</strong> 分頁，該模型另納入診斷年齡、HbA1c 與 eGFR。':''}
        ${(egfr!==null&&egfr<60)||(acr!==null&&acr>=30) ? '<br>依 ESC 2021，中重度 CKD 本身即屬<strong>高／極高心血管風險</strong>；Add-on 的用途是把風險量化得更精準。' : ''}
      </div>
      <table class="regions">
        <thead><tr><th>風險地區</th><th>原始 SCORE2</th><th>＋CKD Add-on</th></tr></thead>
        <tbody>${regionTable(k=>{
            const r = ckdAddon({...p, region:k});
            return [pct(r.base)+'%', r.addon==null?'—':pct(r.addon)+'%'];
          }, c.region)}</tbody>
      </table>`;
  }
  document.querySelectorAll(`#${panelId} select`).forEach(el=>el.addEventListener('change', calc));
  calc();

  document.getElementById(`${pfx}_info`).innerHTML = `
    <h3>這個模型算什麼</h3>
    <ul>
      <li><strong>結果定義</strong>：10 年內第一次發生<strong>致死或非致死</strong>心肌梗塞、中風等心血管事件的機率（已針對競爭死亡風險校正）。</li>
      <li><strong>適用對象</strong>：${kind==='op'?'70–89':'40–69'} 歲、<strong>過去沒有心血管疾病</strong>、也沒有家族性高膽固醇血症的人。</li>
      <li><strong>CKD Add-on 怎麼運作</strong>：先由年齡、性別、血壓、血脂、吸菸、糖尿病推估此人「應該有的 eGFR 與 uACR」，再用實際值與預期值的差距去調整 SCORE2 風險。所以腎功能比同齡族群好的人，風險會被<strong>往下</strong>修正。</li>
      <li><strong>eGFR 未填</strong>時只顯示原始 SCORE2；<strong>只填 eGFR</strong> 時使用 eGFR-only Add-on；兩者都填則用完整版（論文主要模型）。</li>
    </ul>
    <p class="cite">SCORE2${kind==='op'?'-OP':''} working group, <em>Eur Heart J</em> 2021．Matsushita K, et al.（CKD-PC）<em>Eur J Prev Cardiol</em> 2023;30:8-16．</p>`;
}

function initSCORE2DM(){
  buildFields(document.getElementById('g_s2dm'), 's2dm', 'dm');
  const $ = id => document.getElementById(id);
  function calc(){
    const c = readCommon('s2dm');
    const dxage = +$('s2dm_dxage').value;
    const a1c = ngsp2ifcc(+$('s2dm_a1c').value);
    const egfr = +$('s2dm_egfr').value;
    const p = {...c, dxage, a1c, egfr};
    const r = score2dm(p);
    const [cls,label] = clsSCORE2DM(r);
    const dur = c.age - dxage;
    $('s2dm_result').innerHTML = `
      <div class="big-row">
        ${bigBox('10 年心血管風險', pct(r)+'%', cls, `糖尿病病程約 ${dur >= 0 ? dur : 0} 年`, label)}
      </div>
      <div class="interpret">
        <strong>ESC 2023 糖尿病指引分級</strong>：低 &lt; 5%、中度 5–&lt;10%、高 10–&lt;20%、極高 ≥ 20%。
        ${dur<0?'<br><strong>注意</strong>：診斷年齡大於目前年齡，請確認輸入。':''}
      </div>
      <table class="regions">
        <thead><tr><th>風險地區</th><th>10 年心血管風險</th></tr></thead>
        <tbody>${regionTable(k=>[pct(score2dm({...p, region:k}))+'%'], c.region)}</tbody>
      </table>`;
  }
  document.querySelectorAll('#panel-s2dm select').forEach(el=>el.addEventListener('change', calc));
  calc();

  document.getElementById('s2dm_info').innerHTML = `
    <h3>這個模型算什麼</h3>
    <ul>
      <li><strong>適用對象</strong>：40–69 歲、<strong>第 2 型糖尿病</strong>、過去沒有心血管疾病的人。第 1 型糖尿病不適用。</li>
      <li><strong>三個糖尿病專屬變數</strong>：診斷年齡（愈早診斷風險愈高）、HbA1c、eGFR。台灣檢驗報告的 HbA1c 是 %（NGSP），本站會自動換算成論文使用的 mmol/mol（IFCC）。</li>
      <li><strong>已有 ASCVD、嚴重標的器官損傷</strong>（如 eGFR &lt; 45、蛋白尿合併視網膜病變或神經病變）者，本來就屬極高風險，不需再用本模型分級。</li>
      <li>本模型不含 uACR；若想把白蛋白尿一併納入，可參考 SCORE2 + CKD Add-on 分頁（但該版本未針對糖尿病族群專門校正）。</li>
    </ul>
    <p class="cite">SCORE2-Diabetes working group, <em>Eur Heart J</em> 2023;44:2544-56．ESC 2023 Guidelines for the management of cardiovascular disease in patients with diabetes．</p>`;
}

/* ============================================================
   eGFR 換算小工具
   ============================================================ */
function initHelper(){
  const $ = id => document.getElementById(id);
  fill($('h_cr'),  opts(range(0.4,15,0.1), v=>v.toFixed(1)), 1.2);
  fill($('h_age'), opts(range(18,100,1)), 65);
  function go(){
    const e = ckdepi2021(+$('h_cr').value, +$('h_age').value, $('h_sex').value);
    $('h_out').innerHTML = `eGFR ≈ <strong>${e.toFixed(0)}</strong> mL/min/1.73m²（${gBand(e)}）`;
  }
  ['h_cr','h_age','h_sex'].forEach(id=>$(id).addEventListener('change', go));
  go();
}

/* ============================================================
   分頁切換
   ============================================================ */
document.getElementById('tabs').addEventListener('click', e=>{
  const btn = e.target.closest('.tab'); if(!btn) return;
  document.querySelectorAll('.tab').forEach(t=>t.classList.toggle('active', t===btn));
  document.querySelectorAll('.panel').forEach(p=>p.classList.toggle('active', p.id==='panel-'+btn.dataset.panel));
  window.scrollTo({top:0, behavior:'smooth'});
});

initHelper();
initKFRE();
initSCORE2('s2',  'panel-s2',  'std');
initSCORE2('s2op','panel-s2op','op');
initSCORE2DM();
