#include "stdafx.h"
#include "dxB2ShaderResourceStateCache.h"

dx10ShaderResourceStateCache	SRVSManager;

dx10ShaderResourceStateCache::dx10ShaderResourceStateCache()
{
	ResetDeviceState();
}

void dx10ShaderResourceStateCache::ResetDeviceState()
{
	std::memset(m_PSViews, 0, sizeof(m_PSViews));
	std::memset(m_GSViews, 0, sizeof(m_GSViews));
	std::memset(m_VSViews, 0, sizeof(m_VSViews));
#ifdef USE_DX11
    std::memset(m_HSViews, 0, sizeof(m_HSViews));
    std::memset(m_DSViews, 0, sizeof(m_DSViews));
#endif

	m_uiMinPSView = 0xFFFFFFFF;
	m_uiMaxPSView = 0xFFFFFFFF;

	m_uiMinGSView = 0xFFFFFFFF;
	m_uiMaxGSView = 0xFFFFFFFF;

	m_uiMinVSView = 0xFFFFFFFF;
	m_uiMaxVSView = 0xFFFFFFFF;

#ifdef USE_DX11
	m_uiMinHSView = 0xFFFFFFFF;
	m_uiMaxHSView = 0xFFFFFFFF;

	m_uiMinDSView = 0xFFFFFFFF;
	m_uiMaxDSView = 0xFFFFFFFF;
#endif

	m_bUpdatePSViews = false;
	m_bUpdateGSViews = false;
	m_bUpdateVSViews = false;
#ifdef USE_DX11
	m_bUpdateDSViews = false;
	m_bUpdateHSViews = false;
#endif
}

void dx10ShaderResourceStateCache::Apply()
{
	if (m_bUpdatePSViews)
	{
		//HW.pContext->PSSetShaderResources( m_uiMinPSView, m_uiMaxPSView-m_uiMinPSView+1, &m_PSViews[m_uiMinPSView]);
		m_uiMinPSView = 0xFFFFFFFF;
		m_uiMaxPSView = 0xFFFFFFFF;
		m_bUpdatePSViews = false;
	}

	if (m_bUpdateGSViews)
	{
		//HW.pContext->GSSetShaderResources( m_uiMinGSView, m_uiMaxGSView-m_uiMinGSView+1, &m_GSViews[m_uiMinGSView]);
		m_uiMinGSView = 0xFFFFFFFF;
		m_uiMaxGSView = 0xFFFFFFFF;
		m_bUpdateGSViews = false;
	}

	if (m_bUpdateVSViews)
	{
		//HW.pContext->VSSetShaderResources( m_uiMinVSView, m_uiMaxVSView-m_uiMinVSView+1, &m_VSViews[m_uiMinVSView]);
		m_uiMinVSView = 0xFFFFFFFF;
		m_uiMaxVSView = 0xFFFFFFFF;
		m_bUpdateVSViews = false;
	}

#ifdef USE_DX11
	if (m_bUpdateHSViews)
	{
		HW.pContext->HSSetShaderResources( m_uiMinHSView, m_uiMaxHSView-m_uiMinHSView+1, &m_HSViews[m_uiMinHSView]);
		m_uiMinHSView = 0xFFFFFFFF;
		m_uiMaxHSView = 0xFFFFFFFF;
		m_bUpdateHSViews = false;
	}

	if (m_bUpdateDSViews)
	{
		HW.pContext->DSSetShaderResources( m_uiMinDSView, m_uiMaxDSView-m_uiMinDSView+1, &m_DSViews[m_uiMinDSView]);
		m_uiMinDSView = 0xFFFFFFFF;
		m_uiMaxDSView = 0xFFFFFFFF;
		m_bUpdateDSViews = false;
	}

	if (m_bUpdateCSViews)
	{
		HW.pContext->CSSetShaderResources( m_uiMinCSView, m_uiMaxCSView-m_uiMinCSView+1, &m_CSViews[m_uiMinCSView]);
		m_uiMinCSView = 0xFFFFFFFF;
		m_uiMaxCSView = 0xFFFFFFFF;
		m_bUpdateCSViews = false;
	}
#endif
}

void dx10ShaderResourceStateCache::SetPSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxPixelShaderTextures);

	if ( m_PSViews[uiSlot] != pRes)
	{
		m_PSViews[uiSlot] = pRes;
		if (m_bUpdatePSViews)
		{
			m_uiMinPSView = std::min( m_uiMinPSView, uiSlot );
			m_uiMaxPSView = std::max( m_uiMaxPSView, uiSlot );
		}
		else
		{
			m_bUpdatePSViews = true;
			m_uiMinPSView = uiSlot;
			m_uiMaxPSView = uiSlot;
		}
	}
}

void dx10ShaderResourceStateCache::SetGSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxGeometryShaderTextures);

	if ( m_GSViews[uiSlot] != pRes)
	{
		m_GSViews[uiSlot] = pRes;
		if (m_bUpdateGSViews)
		{
			m_uiMinGSView = std::min( m_uiMinGSView, uiSlot );
			m_uiMaxGSView = std::max( m_uiMaxGSView, uiSlot );
		}
		else
		{
			m_bUpdateGSViews = true;
			m_uiMinGSView = uiSlot;
			m_uiMaxGSView = uiSlot;
		}
	}
}

void dx10ShaderResourceStateCache::SetVSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxVertexShaderTextures);

	if ( m_VSViews[uiSlot] != pRes)
	{
		m_VSViews[uiSlot] = pRes;
		if (m_bUpdateVSViews)
		{
			m_uiMinVSView = std::min( m_uiMinVSView, uiSlot );
			m_uiMaxVSView = std::max( m_uiMaxVSView, uiSlot );
		}
		else
		{
			m_bUpdateVSViews = true;
			m_uiMinVSView = uiSlot;
			m_uiMaxVSView = uiSlot;
		}
	}
}

#ifdef USE_DX11
void dx10ShaderResourceStateCache::SetHSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxHullShaderTextures);

	if ( m_HSViews[uiSlot] != pRes)
	{
		m_HSViews[uiSlot] = pRes;
		if (m_bUpdateHSViews)
		{
			m_uiMinHSView = std::min( m_uiMinHSView, uiSlot );
			m_uiMaxHSView = std::max( m_uiMaxHSView, uiSlot );
		}
		else
		{
			m_bUpdateHSViews = true;
			m_uiMinHSView = uiSlot;
			m_uiMaxHSView = uiSlot;
		}
	}
}

void dx10ShaderResourceStateCache::SetDSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxHullShaderTextures);

	if ( m_DSViews[uiSlot] != pRes)
	{
		m_DSViews[uiSlot] = pRes;
		if (m_bUpdateDSViews)
		{
			m_uiMinDSView = std::min( m_uiMinDSView, uiSlot );
			m_uiMaxDSView = std::max( m_uiMaxDSView, uiSlot );
		}
		else
		{
			m_bUpdateDSViews = true;
			m_uiMinDSView = uiSlot;
			m_uiMaxDSView = uiSlot;
		}
	}
}

void dx10ShaderResourceStateCache::SetCSResource( u32 uiSlot, ID3DShaderResourceView	*pRes )
{
	VERIFY(uiSlot<CBackend::mtMaxComputeShaderTextures);

	if ( m_CSViews[uiSlot] != pRes)
	{
		m_CSViews[uiSlot] = pRes;
		if (m_bUpdateCSViews)
		{
			m_uiMinCSView = std::min( m_uiMinCSView, uiSlot );
			m_uiMaxCSView = std::max( m_uiMaxCSView, uiSlot );
		}
		else
		{
			m_bUpdateCSViews = true;
			m_uiMinCSView = uiSlot;
			m_uiMaxCSView = uiSlot;
		}
	}
}
#endif