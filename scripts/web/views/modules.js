import { t } from './shared.js';

const CARDS = [
  {
    track: 'reconstruction',
    icon: '<polygon points="12,2 22,22 2,22" fill="currentColor"/>',
    titleKey: 'home_modules_recon_title',
    descKey: 'home_modules_recon_desc',
    bulletsKey: 'home_modules_recon_feats',
  },
  {
    track: 'calibration',
    icon: '<circle cx="12" cy="12" r="10" fill="currentColor"/>',
    titleKey: 'home_modules_calib_title',
    descKey: 'home_modules_calib_desc',
    bulletsKey: 'home_modules_calib_feats',
  },
  {
    track: 'analysis',
    icon: [
      '<rect x="3" y="3" width="8" height="8" fill="currentColor"/>',
      '<rect x="13" y="3" width="8" height="8" fill="currentColor"/>',
      '<rect x="3" y="13" width="8" height="8" fill="currentColor"/>',
      '<rect x="13" y="13" width="8" height="8" fill="currentColor"/>',
    ].join(''),
    titleKey: 'home_modules_analysis_title',
    descKey: 'home_modules_analysis_desc',
    bulletsKey: 'home_modules_analysis_feats',
  },
  {
    track: 'design',
    icon: '<polygon points="2,2 22,12 2,22" fill="currentColor"/>',
    titleKey: 'home_modules_design_title',
    descKey: 'home_modules_design_desc',
    bulletsKey: 'home_modules_design_feats',
  },
];

export function mountModules(main) {
  main.innerHTML = `
    <section class="page page--modules">
      <header class="page__header">
        <h1 data-i18n="modules_title">ECOMAP Modules</h1>
        <p data-i18n="modules_intro">Four cooperating modules that turn a GEM into a calibrated, design-ready ecGEM. Create or open a project before running them.</p>
      </header>
      <div class="modules-grid">
        ${CARDS.map((c) => `
          <a class="module-card" href="#projects">
            <svg class="module-card__icon module__icon" viewBox="0 0 24 24" aria-hidden="true">${c.icon}</svg>
            <h2 class="module-card__title" data-i18n="${c.titleKey}">${c.titleKey}</h2>
            <p class="module-card__desc" data-i18n="${c.descKey}">${c.descKey}</p>
            <ul class="module-card__bullets" data-bullets-key="${c.bulletsKey}">
              <li data-i18n="${c.bulletsKey}">${c.bulletsKey}</li>
            </ul>
            <span class="module-card__cta" data-i18n="modules_open">Open projects -></span>
          </a>
        `).join('')}
      </div>
    </section>
  `;
  main.querySelectorAll('[data-i18n]').forEach((el) => {
    el.textContent = t(el.getAttribute('data-i18n'), el.textContent);
  });
}
